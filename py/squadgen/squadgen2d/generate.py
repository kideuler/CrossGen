"""
generate.py -- SQuadGen inference: block layouts of planar domain OBJs from a trained SQ-Diffuse (docs/squadgen_2d.md
Sec. 4).

Per domain, Sec. 4's steps as the rest of the package implements them:

1-3. sq_diffuse.domain_condition: the OBJ normalised, N and 4N token points drawn on it, and the frozen Geom-AE's code
     there (N is the checkpoint's own: 2048 for the N2048 runs);
4.   sq_diffuse.sample_latents: SQ-Diffuse's ODE from Taubin-smoothed noise, dopri5 with atol 1e-6 and rtol 1e-3 (the
     release's inference), one seed a sample;
5-7. sq_vae.decode_fields: the frozen SQ-VAE decodes cdf, dcdf, offset and offsetc at the centroids of the extraction
     mesh (roundtrip.extraction_mesh: the OBJ subdivided to about h_rel x its diagonal; 0.002, the release's and the
     round trip's);
8.   QUARRY (quarry.py, in a subprocess: torch's and crossgen's OpenMP runtimes abort in one process once TMOP
     starts): the BlockDecomposition the fields describe, meshed by BlockQuadMesh and smoothed by TMOP.

SQ-Diffuse is meant to produce a good layout of the domain, not a particular one, and not every sample extracts: a chart
QUARRY cannot make a quad, or a wall with a gap, leaves blocks that do not cover the domain. So each domain gets
--samples seeds and keeps the best of them: blocks that cover the domain first, then a mesh with no inverted quad, then
the highest alignment quality (BlockDecomposition.alignment_quality(), the metric ORACLE's selector ranks methods by).
A sample is `valid` when it passes the first two.

Given the CDF data the model was trained on (--data; by default the checkpoint's data directory, when it is here), each
domain's layout there is drawn beside the samples, and its block qualities (the manifest's, which cdf-data computes with
the same BlockDecomposition methods QUARRY's blocks are scored by) are reported beside theirs. A domain is `held out`
when its family is one the SQ-VAE's Geom-AE held out (recorded in the SQ-VAE checkpoint), so that no network in the
chain trained on it; the others were training domains, whose data layout SQ-Diffuse has seen.

Device: the first of CUDA, MPS and the CPU, as the trainers pick. Training at N = 2048 is a CUDA job, but sampling one
latent at a time is not: on an 8 GB M2 a sample takes ~15 s on MPS (~42 s on the CPU) and matches the CPU's. A latent
that comes back non-finite is drawn again on the CPU, and the rest of the run stays there.

Run
---
    python -m squadgen2d generate data/meshes/singlemat/geom002.obj                    # -> py/squadgen/runs/generate
    python -m squadgen2d generate data/meshes/singlemat/*.obj --data ~/Documents/2d_CDF_data --samples 4
    python -m squadgen2d generate my_domain.obj --model py/squadgen/runs/<run>/best.pt --out /tmp/gen

Outputs in --out:
* samples.csv: one row per sample (seed, ODE function evaluations, seconds, whether the blocks cover the domain, blocks,
  charts that made no block, inverted quads and min scaled Jacobian of the smoothed mesh, the four block qualities,
  valid, kept);
* domains.csv: one row per domain (held out or training, valid samples, the kept sample and its qualities, and the data
  layout's method, blocks and qualities);
* <model>/seed<k>/: latent.npy (N x C, denormalised), sample.json (seed, evaluations, seconds, QUARRY's exit code)
  and quarry_output.npz (quarry.py: chart of every extraction-mesh triangle, the blocks, and the mesh); with --keep,
  quarry_input.npz too (the mesh and the decoded fields);
* <model>/quad.ply: the kept sample's smoothed quad mesh (the best one that covers, valid or not: domains.csv says);
* generate.png (generate_<page>.png past --rows domains): a row per domain -- the data layout with its chart sides,
  every sample's decoded CDF with QUARRY's block edges, and the kept sample's mesh;
* log.txt.

--replot redoes the CSVs and figures from those files (the latents decoded again, nothing sampled or extracted), so a
change to the figure does not cost the sampling again.
"""

import argparse
import csv
import json
import os
import shutil
import time

import numpy as np
import torch

from . import training as tr
from .domains import family
from .layout import Layout, chart_sides, layout_fields, read_manifest
from .paths import REPO, RUNS, from_repo, to_repo
from .quarry import DEFAULT_SWEEPS, QUALITIES
from .roundtrip import extraction_mesh, locate, run_quarry, write_ply, write_quarry_input
from .sq_diffuse import cond_to, domain_condition, load_sq_diffuse, sample_latents
from .sq_vae import FIELD_RAMP, decode_fields

DEFAULT_MODEL = os.path.join(REPO, "best_models", "best.pt")
DEFAULT_OUT = os.path.join(RUNS, "generate")
SAMPLE_KEYS = ("seed", "nfe", "seconds", "covers", "blocks", "failed", "inverted", "min_sj") + QUALITIES + (
    "valid", "kept")


class Generator:
    """The checkpoint's three networks on one device, moved to the CPU for good if a latent comes back non-finite."""

    def __init__(self, path, dev, log, vae_path=None):
        self.path, self.dev, self.log = path, dev, log
        self.sit, self.vae, self.geom, self.cfg, self.vcfg, self.gcfg = load_sq_diffuse(path, dev, vae_path)

    def to_cpu(self):
        self.log(f"  a non-finite latent on {self.dev}: drawing it again on the CPU, and the rest of the run there")
        self.dev = "cpu"
        for m in (self.sit, self.vae, self.geom):
            m.to("cpu")

    def condition(self, obj, seed):
        """domain_condition on the CPU-held copy: (cond, Domain, boundary segments)."""
        cond, d, bnd = domain_condition(self.geom, self.gcfg, obj, self.dev, seed=seed)
        return {k: v.cpu() for k, v in cond.items()}, d, bnd

    @torch.no_grad()
    def latent(self, cond, bnd, seed, args):
        z, nfe = sample_latents(self.sit, cond_to(cond, self.dev), [seed], args.noise, [bnd], args.method, args.steps,
                                args.atol, args.rtol)
        if not torch.isfinite(z).all():
            if self.dev == "cpu":
                raise RuntimeError(f"seed {seed}: the latent is not finite on the CPU either")
            self.to_cpu()
            return self.latent(cond, bnd, seed, args)
        return z, nfe

    @torch.no_grad()
    def fields(self, z, P):
        return decode_fields(self.vae, z, P)


def read_output(d):
    """What a sample's quarry_output.npz says, as a flat dict, with the block edges as polylines."""
    out = {"covers": False, "blocks": 0, "failed": -1, "inverted": -1, "min_sj": float("nan"), "edges": [],
           "mesh": None, **{q: float("nan") for q in QUALITIES}}
    path = os.path.join(d, "quarry_output.npz")
    if not os.path.exists(path):
        return out
    with np.load(path) as z:
        rep = json.loads(str(z["report"]))
        pts, start = z["edge_points"], z["edge_start"]
        out.update(covers=bool(z["covers"]), blocks=len(z["blocks"]), failed=int(z["failed"]),
                   edges=[pts[a:b] for a, b in zip(start[:-1], start[1:]) if b - a > 1])
        if "mesh_quads" in z.files:
            out.update(inverted=int(z["inverted"]), min_sj=float(z["min_sj"]),
                       mesh=(z["mesh_vertices"].copy(), z["mesh_quads"].copy()))
    out.update({q: float(rep[q]) for q in QUALITIES if q in rep})
    return out


def rank(s):
    """The order samples are kept in: blocks that cover the domain, an untangled mesh, then alignment quality."""
    aq = s["alignment_quality"]
    return (s["covers"], s["covers"] and s["inverted"] == 0, aq if np.isfinite(aq) else -1.0)


def reference(rows, obj, V, T):
    """The domain's layout in the CDF data, on the extraction mesh: (row, cdf per triangle, chart sides in OBJ units),
    or None. A domain with several layouts (one per method) shows the first."""
    rs = rows.get(os.path.splitext(os.path.basename(obj))[0]) if rows else None
    if not rs:
        return None
    r = rs[0]
    L = Layout(r["npz"])
    tri, bary, _ = locate(L.V, L.T, (V[T].mean(1) - L.centre) * L.scale)
    return {"row": r, "cdf": layout_fields(L, L.V, tri, bary)["cdf"][:, 0],
            "sides": chart_sides(r["npz"], L.V) / L.scale + L.centre}


def failure(s):
    """Why a sample's blocks do not cover the domain, in a few words."""
    if s["rc"] != 0:
        return f"QUARRY failed (rc {s['rc']})"
    if s["failed"] > 0:
        return f"{s['failed']} charts/sides unmatched"
    return "blocks not a valid decomposition"


def generate_domain(gen, obj, args, val_fams, rows, log):
    """Sample one domain (or, with --replot, re-read its samples from --out) and keep the best."""
    name = os.path.splitext(os.path.basename(obj))[0]
    held = family(obj) in val_fams if val_fams is not None else None
    t0 = time.time()
    cond, d, bnd = gen.condition(obj, args.seed)
    V, T, levels = extraction_mesh(obj, args.h_rel)
    P = ((V[T].mean(1) - d.centre) * d.scale).astype(np.float32)
    dom_dir = os.path.join(args.out, name)
    samples = []
    for k in range(args.samples):
        seed = args.seed + k
        sd = os.path.join(dom_dir, f"seed{seed}")
        t = time.time()
        if args.replot:
            if not os.path.exists(os.path.join(sd, "sample.json")):
                raise SystemExit(f"--replot: no sample of {name} seed {seed} in {to_repo(dom_dir)}")
            z = torch.as_tensor(np.load(os.path.join(sd, "latent.npy")))[None].to(gen.dev)
            f = gen.fields(z, P)
            with open(os.path.join(sd, "sample.json")) as fh:
                meta = json.load(fh)
        else:
            z, nfe = gen.latent(cond, bnd, seed, args)
            f = gen.fields(z, P)
            shutil.rmtree(sd, ignore_errors=True)
            write_quarry_input(sd, V, T, f, d.scale, args.sweeps)
            np.save(os.path.join(sd, "latent.npy"), z[0].cpu().numpy())
            rc = run_quarry(None, sd)["rc"]
            meta = {"seed": seed, "nfe": nfe, "seconds": round(time.time() - t, 1), "rc": rc}
            with open(os.path.join(sd, "sample.json"), "w") as fh:
                json.dump(meta, fh)
            if not args.keep:
                os.remove(os.path.join(sd, "quarry_input.npz"))
        s = read_output(sd)
        s.update(meta, cdf=f["cdf"][:, 0])
        s["valid"] = bool(s["covers"] and s["inverted"] == 0)
        samples.append(s)
        log(f"  {name} seed {seed}: " + (f"{s['blocks']} blocks, " + ("valid" if s["valid"] else
                                                                      f"{s['inverted']} quads inverted")
                                         + f", alignment {s['alignment_quality']:.3f}" if s["covers"] else
                                         f"not extracted: {failure(s)}")
            + f"  ({s['nfe']} evaluations, {s['seconds']} s)")
    best = max(samples, key=rank)
    kept = best if best["covers"] else None
    for s in samples:
        s["kept"] = s is kept
    if kept is not None and kept["mesh"] is not None and not args.replot:
        write_ply(os.path.join(dom_dir, "quad.ply"), *kept["mesh"])
    ref = reference(rows, obj, V, T)
    out = {"name": name, "obj": obj, "held": held, "V": V, "T": T, "samples": samples, "kept": kept, "ref": ref,
           "seconds": time.time() - t0}
    log(f"{name} ({'held out' if held else 'training' if held is not None else 'split unknown'}): "
        f"{sum(s['valid'] for s in samples)}/{len(samples)} valid; "
        + (f"kept seed {kept['seed']}, {kept['blocks']} blocks, alignment {kept['alignment_quality']:.3f}"
           + ("" if kept["valid"] else f" (not valid: {kept['inverted']} quads inverted)") if kept else "none extracted")
        + (f"; data: {ref['row']['method']}, {ref['row']['num_blocks']} blocks, alignment "
           f"{float(ref['row']['alignment_quality']):.3f}" if ref else "")
        + f"  ({out['seconds']:.0f} s)")
    return out


def plot(results, out, title):
    """A row per domain: the data layout (when there is one) with its chart sides, every sample's decoded CDF with
    QUARRY's block edges, and the kept sample's smoothed mesh with its block edges. Every panel is in OBJ coordinates;
    a sample's state is written in its title, not only shown by its lines."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap

    ink, muted, edge = "#1a1a19", "#52514e", "#c9c8c4"
    cmap = LinearSegmentedColormap.from_list("field", FIELD_RAMP)
    has_ref = any(r["ref"] is not None for r in results)
    k = max(len(r["samples"]) for r in results)
    cols = int(has_ref) + k + 1
    aspect = [float(np.clip(np.ptp(r["V"][:, 1]) / np.ptp(r["V"][:, 0]), 0.45, 1.3)) for r in results]
    fig, axes = plt.subplots(len(results), cols, figsize=(2.5 * cols, sum(2.25 * a + 0.75 for a in aspect) + 0.6),
                             squeeze=False, gridspec_kw={"height_ratios": [a + 0.33 for a in aspect]})
    im = None

    def field(ax, V, T, vals):
        return ax.tripcolor(V[:, 0], V[:, 1], T, facecolors=np.clip(vals, 0, 1), cmap=cmap, vmin=0, vmax=1,
                            shading="flat", edgecolors="none", rasterized=True)

    def lines(ax, segs, width, colour=ink):
        ax.add_collection(LineCollection(segs, colors=colour, linewidths=width))

    for row, r in zip(axes, results):
        V, T = r["V"], r["T"]
        split = {True: "held out", False: "training", None: ""}[r["held"]]
        j = 0
        if has_ref:
            ax = row[0]
            if r["ref"] is not None:
                rr = r["ref"]["row"]
                im = field(ax, V, T, r["ref"]["cdf"])
                lines(ax, r["ref"]["sides"], 0.6)
                ax.set_title(f"{r['name']} ({split})\ndata: {rr['method']}, {rr['num_blocks']} blocks, "
                             f"align {float(rr['alignment_quality']):.2f}", fontsize=7.5, color=muted)
            else:
                ax.set_title(f"{r['name']} ({split})\nnot in the data", fontsize=7.5, color=muted)
            j = 1
        for s in r["samples"]:
            ax = row[j]
            im = field(ax, V, T, s["cdf"])
            if s["edges"]:
                lines(ax, s["edges"], 0.8)
            if s["covers"]:
                state = ("valid" if s["inverted"] == 0 else f"{s['inverted']} inverted") + \
                        f", align {s['alignment_quality']:.2f}"
                head = f"seed {s['seed']}: {s['blocks']} blocks" + (
                    (" (kept)" if s["valid"] else " (best)") if s["kept"] else "")
            else:
                head, state = f"seed {s['seed']}: not extracted", failure(s)
            if j == 0:
                head = f"{r['name']} ({split})  {head}"
            ax.set_title(f"{head}\n{state}", fontsize=7.5, color=ink if s["kept"] else muted)
            j += 1
        for ax in row[j:-1]:
            ax.axis("off")
        ax = row[-1]
        kept = r["kept"]
        if kept is not None and kept["mesh"] is not None:
            QV, Q = kept["mesh"]
            quad_edges = QV[np.concatenate([Q[:, [0, 1]], Q[:, [1, 2]], Q[:, [2, 3]], Q[:, [3, 0]]])]
            lines(ax, quad_edges, 0.25, "#8a8984")
            lines(ax, kept["edges"], 0.9)
            if kept["valid"]:
                ax.set_title(f"kept: seed {kept['seed']}, {len(Q)} quads\nmin scaled Jacobian {kept['min_sj']:.2f}",
                             fontsize=7.5, color=ink)
            else:
                ax.set_title(f"best, NOT valid: seed {kept['seed']}\n{kept['inverted']} of {len(Q)} quads inverted",
                             fontsize=7.5, color=muted)
        else:
            ax.set_title("no sample made blocks\nthat cover the domain", fontsize=7.5, color=muted)
        ax.set_xlim(V[:, 0].min(), V[:, 0].max())
        ax.set_ylim(V[:, 1].min(), V[:, 1].max())
        for ax in row:
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color(edge)
    fig.suptitle(title, y=0.997, fontsize=9, color=ink)
    fig.tight_layout(rect=(0, 0, 0.92, 0.975))
    if im is not None:
        cax = fig.add_axes((0.935, 0.3, 0.01, 0.4))
        cb = fig.colorbar(im, cax=cax, ticks=[0, 0.25, 0.5, 0.75, 1])
        cb.set_label("CDF (data: true; samples: decoded)", fontsize=7.5)
        cb.ax.tick_params(labelsize=7)
    fig.savefig(out, dpi=110, facecolor="#fcfcfb")
    plt.close(fig)


def write_csvs(results, out):
    with open(os.path.join(out, "samples.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("model", "split") + SAMPLE_KEYS)
        for r in results:
            for s in r["samples"]:
                w.writerow([r["name"], {True: "held_out", False: "training", None: ""}[r["held"]]]
                           + [round(s[k], 4) if isinstance(s[k], float) else s[k] for k in SAMPLE_KEYS])
    keys = ("blocks",) + QUALITIES
    with open(os.path.join(out, "domains.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "split", "samples", "valid", "kept_seed"] + [f"kept_{k}" for k in keys]
                   + ["data_method"] + [f"data_{k}" for k in keys])
        for r in results:
            kept, ref = r["kept"], r["ref"]
            w.writerow([r["name"], {True: "held_out", False: "training", None: ""}[r["held"]], len(r["samples"]),
                        sum(s["valid"] for s in r["samples"]), kept["seed"] if kept else ""]
                       + [(round(kept[k], 4) if isinstance(kept[k], float) else kept[k]) if kept else "" for k in keys]
                       + ([ref["row"]["method"], ref["row"]["num_blocks"]]
                          + [round(float(ref["row"][q]), 4) for q in QUALITIES] if ref else [""] * (len(keys) + 1)))


def summarise(results, log):
    """Valid samples and domains, and the kept layouts' qualities against the data's, overall and per split."""
    groups = [("all", results)] + [(n, [r for r in results if r["held"] is h]) for n, h in (("held out", True),
                                                                                          ("training", False))]
    for name, rs in groups:
        if not rs:
            continue
        ss = [s for r in rs for s in r["samples"]]
        both = [r for r in rs if r["kept"] is not None and r["kept"]["valid"] and r["ref"] is not None]
        line = (f"{name}: {len(rs)} domains; samples valid {sum(s['valid'] for s in ss)}/{len(ss)}, cover "
                f"{sum(s['covers'] for s in ss)}/{len(ss)}; domains with a valid sample "
                f"{sum(any(s['valid'] for s in r['samples']) for r in rs)}/{len(rs)}")
        if both:
            line += f"; kept (valid) vs data on those {len(both)}:"
            for q in QUALITIES:
                line += (f" {q} {np.mean([r['kept'][q] for r in both]):.3f}/"
                         f"{np.mean([float(r['ref']['row'][q]) for r in both]):.3f}")
            line += (f", blocks {np.mean([r['kept']['blocks'] for r in both]):.1f}/"
                     f"{np.mean([int(r['ref']['row']['num_blocks']) for r in both]):.1f}")
        log(line)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="squadgen2d generate", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("objs", nargs="+", help="domain OBJ files (planar triangle meshes)")
    ap.add_argument("--model", default=DEFAULT_MODEL, help="the SQ-Diffuse checkpoint (default best_models/best.pt)")
    ap.add_argument("--vae", help="its SQ-VAE, when not at the path the checkpoint records")
    ap.add_argument("--data", help="the CDF data it was trained on, for the reference layouts (default: the "
                    "checkpoint's data directory, when it is here)")
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--samples", type=int, default=4, help="seeds per domain; the best valid one is kept")
    ap.add_argument("--seed", type=int, default=0, help="first seed (sample k uses seed + k)")
    ap.add_argument("--noise", default="taubin", choices=("taubin", "plain"),
                    help="initial noise: Taubin-smoothed (the release's inference) or white (as trained)")
    ap.add_argument("--method", default="dopri5", choices=("dopri5", "heun", "euler"))
    ap.add_argument("--steps", type=int, default=12, help="steps of a fixed-step --method")
    ap.add_argument("--atol", type=float, default=1e-6)
    ap.add_argument("--rtol", type=float, default=1e-3)
    ap.add_argument("--h_rel", type=float, default=0.002, help="extraction-mesh edge / bbox diagonal (release: 0.002)")
    ap.add_argument("--sweeps", type=int, default=DEFAULT_SWEEPS, help="TMOP sweeps on the block mesh")
    ap.add_argument("--rows", type=int, default=6, help="domains per figure page (0: no figure)")
    ap.add_argument("--keep", action="store_true", help="keep each sample's quarry_input.npz (mesh and fields)")
    ap.add_argument("--replot", action="store_true",
                    help="redo the CSVs and figures of an earlier run in --out from its saved samples (latent.npy, "
                    "quarry_output.npz), without sampling or extracting again; same OBJs, --samples and --seed")
    ap.add_argument("--device", default="auto", help="cuda, cuda:N, mps or cpu (default: the first there is)")
    args = ap.parse_args(argv)
    if args.samples < 1:
        raise SystemExit("--samples must be at least 1")
    objs = [os.path.abspath(o) for o in args.objs]
    missing = [o for o in objs if not os.path.exists(o)]
    if missing:
        raise SystemExit(f"no such OBJ: {missing[0]}" + (f" (and {len(missing) - 1} more)" if len(missing) > 1 else ""))
    os.makedirs(args.out, exist_ok=True)
    logf = open(os.path.join(args.out, "log.txt"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    dev = tr.pick_device(args.device)
    gen = Generator(args.model, dev, log, args.vae)
    cfg = gen.cfg
    data = args.data or (from_repo(cfg["data"]) if cfg.get("data") else None)
    rows = None
    if data and os.path.exists(os.path.join(data, "manifest.csv")):
        rows = {}
        for r in read_manifest(data):
            rows.setdefault(r["model"], []).append(r)
    elif args.data:
        raise SystemExit(f"no manifest.csv in --data {args.data}")
    vf = gen.vcfg.get("val_families")
    val_fams = set(vf) if vf is not None else None
    ck_step = torch.load(args.model, map_location="cpu", weights_only=False)["step"]
    log(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  {'replot of ' if args.replot else ''}{to_repo(args.model)} "
        f"(step {ck_step}, N = {cfg.get('tokens')}), "
        f"SQ-VAE {cfg['vae']}, device {dev}; {len(objs)} domains x {args.samples} samples, {args.noise} noise, "
        f"{args.method}, h_rel {args.h_rel}; reference layouts: {to_repo(data) if rows else 'none'}")
    results = []
    for obj in objs:
        results.append(generate_domain(gen, obj, args, val_fams, rows, log))
    write_csvs(results, args.out)
    summarise(results, log)
    if args.rows:
        pages = [results[i:i + args.rows] for i in range(0, len(results), args.rows)]
        title = (f"SQ-Diffuse {to_repo(args.model)} (step {ck_step}, N = {cfg.get('tokens')}): {args.samples} samples "
                 f"per domain through QUARRY; lines: chart sides (data) / block edges (samples)")
        for n, page in enumerate(pages, 1):
            path = os.path.join(args.out, "generate.png" if len(pages) == 1 else f"generate_{n}.png")
            plot(page, path, title)
            log(f"figure -> {to_repo(path)}")
    logf.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
