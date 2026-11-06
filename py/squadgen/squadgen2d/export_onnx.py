"""
export_onnx.py -- the trained SQuadGen networks as ONNX graphs, for inference without Python (src/SQUADGEN).

Sec. 4 inference runs three networks: the Geom-AE's encoder (the geometry code at the tokens), SQ-Diffuse's SiT (the
velocity its ODE integrates, about 110 times a sample) and the SQ-VAE's decoder (the fields at any point). They are
exported as five graphs, cut where the computation has a seam, so that nothing is recomputed that has not changed:

    geom_code.onnx      xyz (1, N, 2), sdf (1, N, 1), xyz_kv (1, 4N, 2), sdf_kv (1, 4N, 1) -> code (1, N, code_dim):
                        the Geom-AE's return_last_features_and_norm of its posterior mean, as domain_condition takes
                        it. Once per domain.
    sit_condition.onnx  code, xyz, sdf -> pe (depth, 1, N, hidden): what each SiT block adds to the tokens before it
                        (SiT.cond_embed). It depends on the geometry alone, so a sampler computes it once per domain,
                        as sample_latents does, and not once per function evaluation.
    sit_velocity.onnx   x (1, N, C), t (1,), pe -> v (1, N, C): one function evaluation of the ODE, on normalised
                        latents (target_mu / target_std are in the metadata; the latent the decoder reads is
                        x * target_std + target_mu).
    vae_tokens.onnx     z (1, N, C) -> tokens (1, N, vae_dim): the SQ-VAE decoder's proj and self-attention
                        (return_last_features). Once per latent, where decode_fields repeats it for every chunk of
                        query points.
    vae_query.onnx      tokens, xy (1, P, 2) -> cdf (1, P, 1), dcdf (1, P, 1), offset (1, P, 2), offsetc (1, P, 2):
                        the four fields QUARRY reads (query_last_features, offset divided by 7 as there), at any P.
                        The other five output heads are left out: nothing downstream reads them.

Every point is in the domain's normalised frame (geometry.normalisation: area centroid at the origin, farthest
boundary vertex at radius 0.999) and every field comes out in it: the offsets are frame vectors, which a caller
divides by the frame's scale for the model's units, as write_quarry_input does. N, 4N and P are dynamic axes; the
batch is 1. Everything is float32, as the networks were trained: the SQ-VAE's decoder tolerates very little latent
error (memory squadgen-sq-diffuse-2d), so nothing is quantised or run in half precision.

What a C++ reader needs besides the weights is in each file's metadata, under squadgen.* keys: N (`tokens`), the
latent width, the SiT's depth and width, the code and decoder widths, tau (the truncated SDF's), target_mu and
target_std, and the checkpoint and step it was exported from. All five carry the same `export` id (the SiT's and the
SQ-VAE's weight fingerprints), so a reader can refuse a directory that mixes two exports. squadgen.json repeats it for
people.

Two changes to the networks, neither of them to what they compute:
* nn.RMSNorm (the SiT's QK-norm) has no rule in torch's TorchScript exporter, which is the one used here because the
  torch.export one needs onnxscript. It is replaced, in a copy, by its definition: x / sqrt(mean(x^2) + eps) * weight,
  eps float32's machine epsilon as nn.RMSNorm takes it when given none.
* The graphs are the modules' own forward passes, wrapped to take tensors rather than the dicts the trainers pass.

The reference case. The export also runs Sec. 4 in torch on one domain (--reference, by default singlemat geom005:
eight blocks, valid on every seed of the step-285k model) and writes reference/: the token sets and their sdf, the
geometry code, the white noise of the seed and its Taubin-smoothed form, the velocity there at t = 0, the latent
dopri5 reaches and how many evaluations it took, the decoder tokens, the extraction mesh and the decoded fields at its
centroids, and what QUARRY made of them. `TestSQUADGEN --check DIR/reference` runs the C++ pipeline from the same
inputs and compares every stage with these, which is the test that the C++ port and the graphs are the networks.

When onnxruntime is importable, each graph is also run here on the reference inputs and compared with torch. It is
not a dependency: the C++ check covers the same ground with the runtime the C++ links.

Run
---
    python -m squadgen2d export-onnx                                     # best_models/best.pt -> best_models/onnx
    python -m squadgen2d export-onnx --model py/squadgen/runs/<run>/best.pt --out /tmp/sqd_onnx
    python -m squadgen2d export-onnx --reference data/meshes/singlemat/geom002.obj --seed 3
"""

import argparse
import copy
import hashlib
import json
import os
import shutil
import time
import warnings

import numpy as np
import torch
import torch.nn as nn

from . import geometry as geo
from .domains import Domain
from .geom_ae import make_sample as geom_sample
from .generate import DEFAULT_MODEL, read_output
from .paths import REPO, to_repo
from .quarry import DEFAULT_SWEEPS
from .roundtrip import extraction_mesh, run_quarry, write_quarry_input
from .sq_diffuse import TAUBIN, dopri5, domain_condition, fingerprint, load_sq_diffuse, sample_latents, taubin_smooth
from .sq_vae import decode_fields, geom_code

DEFAULT_OUT = os.path.join(REPO, "best_models", "onnx")
DEFAULT_REFERENCE = os.path.join(REPO, "data", "meshes", "singlemat", "geom005.obj")
GRAPHS = ("geom_code", "sit_condition", "sit_velocity", "vae_tokens", "vae_query")
QUERY_KEYS = ("cdf", "dcdf", "offset", "offsetc")
OPSET = 17      # LayerNormalization is an operator from 17 on, which ONNX Runtime runs fused


# =============================================================================
# 1. THE GRAPHS
# =============================================================================


class ExplicitRMSNorm(nn.Module):
    """nn.RMSNorm over the last dimension, written out: the TorchScript exporter has no rule for aten::rms_norm."""

    def __init__(self, norm):
        super().__init__()
        self.weight = norm.weight
        self.eps = norm.eps

    def forward(self, x):
        eps = torch.finfo(x.dtype).eps if self.eps is None else self.eps
        return x * torch.rsqrt(x.pow(2).mean(-1, keepdim=True) + eps) * self.weight


class GeomCode(nn.Module):
    """The Geom-AE's code at the tokens: geom_code() with the token sets as tensors (KLAutoEncoder.encode at its
    posterior mean, then return_last_features_and_norm)."""

    def __init__(self, geom):
        super().__init__()
        self.geom = geom

    def forward(self, xyz, sdf, xyz_kv, sdf_kv):
        g = self.geom
        x = g.point_embed(xyz, {"sdf": sdf})
        pc = g.point_embed(xyz_kv, {"sdf": sdf_kv})
        cross_attn, cross_ff = g.cross_attend_blocks
        x = cross_attn(x, context=pc) + x
        x = cross_ff(x) + x
        return g.return_last_features_and_norm(g.mean_fc(x))


class SiTCondition(nn.Module):
    """SiT.cond_embed, the per-block conditioning, stacked: (depth, 1, N, hidden)."""

    def __init__(self, sit):
        super().__init__()
        self.sit = sit

    def forward(self, code, xyz, sdf):
        return torch.stack(self.sit.cond_embed({"code": code, "xyz": xyz, "sdf": sdf}))


class SiTVelocity(nn.Module):
    """SiT.forward with the conditioning given (pe), as sample_latents calls it."""

    def __init__(self, sit):
        super().__init__()
        self.sit = sit

    def forward(self, x, t, pe):
        s = self.sit
        x = s.input_proj(x)
        c = s.t_embedder(t)
        for i, block in enumerate(s.blocks):
            x = block(x + pe[i], c)
        return s.final_layer(x, c)


class VaeTokens(nn.Module):
    def __init__(self, vae):
        super().__init__()
        self.vae = vae

    def forward(self, z):
        return self.vae.return_last_features(z)


class VaeQuery(nn.Module):
    """SQVAE.query_last_features for the four keys QUARRY reads."""

    def __init__(self, vae):
        super().__init__()
        self.vae = vae

    def forward(self, tokens, xy):
        v = self.vae
        h = v.decoder_cross_attn(v.point_embed.query(xy), context=tokens)
        out = v.to_outputs
        return out["cdf"](h), out["dcdf"](h), out["offset"](h) / 7, out["offsetc"](h)


def graph_specs(sit, vae, geom):
    """name -> (module, input names, output names, dynamic axes). The SiT is a copy with its RMSNorms written out."""
    sit = copy.deepcopy(sit)
    for block in sit.blocks:
        block.attn.q_norm = ExplicitRMSNorm(block.attn.q_norm)
        block.attn.k_norm = ExplicitRMSNorm(block.attn.k_norm)
    n, m, p = {1: "n"}, {1: "m"}, {1: "p"}
    return {
        "geom_code": (GeomCode(geom), ["xyz", "sdf", "xyz_kv", "sdf_kv"], ["code"],
                      {"xyz": n, "sdf": n, "xyz_kv": m, "sdf_kv": m, "code": n}),
        "sit_condition": (SiTCondition(sit), ["code", "xyz", "sdf"], ["pe"],
                          {"code": n, "xyz": n, "sdf": n, "pe": {2: "n"}}),
        "sit_velocity": (SiTVelocity(sit), ["x", "t", "pe"], ["v"], {"x": n, "pe": {2: "n"}, "v": n}),
        "vae_tokens": (VaeTokens(vae), ["z"], ["tokens"], {"z": n, "tokens": n}),
        "vae_query": (VaeQuery(vae), ["tokens", "xy"], list(QUERY_KEYS),
                      {"tokens": n, "xy": p, **{k: p for k in QUERY_KEYS}}),
    }


def example_inputs(sit, vae, geom, n):
    """Inputs of the right shapes to trace each graph with (a small P: the axis is dynamic)."""
    g = torch.Generator().manual_seed(0)
    xyz = torch.rand((1, n, 2), generator=g) - 0.5
    sdf = -torch.rand((1, n, 1), generator=g)
    xyz_kv = torch.rand((1, 4 * n, 2), generator=g) - 0.5
    sdf_kv = -torch.rand((1, 4 * n, 1), generator=g)
    with torch.no_grad():
        code = GeomCode(geom)(xyz, sdf, xyz_kv, sdf_kv)
        pe = SiTCondition(sit)(code, xyz, sdf)
        z = torch.randn((1, n, sit.channels), generator=g)
        tokens = vae.return_last_features(sit.denormalize(z))
    return {"geom_code": (xyz, sdf, xyz_kv, sdf_kv), "sit_condition": (code, xyz, sdf),
            "sit_velocity": (z, torch.tensor([0.5]), pe), "vae_tokens": (sit.denormalize(z),),
            "vae_query": (tokens, torch.rand((1, 1000, 2), generator=g) - 0.5)}


def metadata(sit, vae, cfg, gcfg, model_path, step, vae_fp):
    """The squadgen.* metadata every graph carries (see the module docstring)."""
    sit_fp = fingerprint(sit)
    fmt = lambda xs: ",".join(format(float(x), ".9g") for x in xs)   # float32 round-trips at 9 digits
    return {
        "squadgen.export": hashlib.sha1(f"{sit_fp}:{vae_fp}".encode()).hexdigest()[:16],
        "squadgen.tokens": str(int(cfg["tokens"])),
        "squadgen.latent_dim": str(int(sit.channels)),
        "squadgen.depth": str(len(sit.blocks)),
        "squadgen.hidden": str(int(sit.input_proj.out_features)),
        "squadgen.code_dim": str(int(cfg["in_channels_cond"])),
        "squadgen.vae_dim": str(int(vae.proj.out_features)),
        "squadgen.tau": format(float(gcfg["tau"]), ".17g"),
        "squadgen.target_mu": fmt(sit.target_mu.tolist()),
        "squadgen.target_std": fmt(sit.target_std.tolist()),
        "squadgen.model": to_repo(model_path),
        "squadgen.step": str(int(step)),
        "squadgen.vae": cfg["vae"],
        "squadgen.vae_fingerprint": vae_fp,
        "squadgen.sit_fingerprint": sit_fp,
        "squadgen.frame": "area centroid at the origin, farthest boundary vertex at radius 0.999",
        "squadgen.query_outputs": ",".join(QUERY_KEYS),
    }


def export_graphs(sit, vae, geom, meta, out, opset, log):
    """Every graph to out/<name>.onnx with the metadata; returns {name: seconds}."""
    import onnx

    specs = graph_specs(sit, vae, geom)
    examples = example_inputs(sit, vae, geom, int(meta["squadgen.tokens"]))
    times = {}
    for name in GRAPHS:
        module, ins, outs, dyn = specs[name]
        path = os.path.join(out, f"{name}.onnx")
        t0 = time.time()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            torch.onnx.export(module.eval(), examples[name], path, dynamo=False, opset_version=opset,
                              input_names=ins, output_names=outs, dynamic_axes=dyn, do_constant_folding=True)
        m = onnx.load(path)
        for k, v in dict(meta, **{"squadgen.graph": name}).items():
            entry = m.metadata_props.add()
            entry.key, entry.value = k, v
        m.doc_string = f"SQuadGen 2D {name} (python -m squadgen2d export-onnx; src/SQUADGEN reads it)"
        onnx.save(m, path)
        times[name] = time.time() - t0
        log(f"  {name}.onnx: {os.path.getsize(path) / 1e6:.1f} MB ({times[name]:.1f} s)")
    return times


# =============================================================================
# 2. THE REFERENCE CASE
# =============================================================================


def condition_inputs(geom, gcfg, obj, seed):
    """domain_condition's token sets and their sdf, drawn as it draws them when no cache holds the domain: (Domain,
    {xyz, sdf, xyz_kv, sdf_kv} as (1, k, c) float32 tensors)."""
    d = Domain(obj)
    n = geom.num_latents
    rng = np.random.default_rng(seed)
    d.fps = {}
    for m in (n, 4 * n):
        tri, bary = geo.token_sets(d.V, d.F, d.cdf, m, 1, rng)
        d.fps[m] = (tri, bary.astype(np.float32))
    s = geom_sample(d, dict(gcfg, n_tokens=n, n_query=1, n_near=1, n_global=1), np.random.default_rng(seed), False)
    t = {k: torch.from_numpy(s[f"{k}_fps_{m}"])[None] for k, m in (("xyz", n), ("sdf", n))}
    t.update({f"{k}_kv": torch.from_numpy(s[f"{k}_fps_{4 * n}"])[None] for k in ("xyz", "sdf")})
    return d, t


@torch.no_grad()
def reference_case(sit, vae, geom, gcfg, obj, seed, h_rel, atol, rtol, out, log):
    """Sec. 4 in torch on one domain, every intermediate saved to out/ (see the module docstring). Returns the
    numbers case.txt holds, and the arrays, for the onnxruntime check."""
    os.makedirs(out, exist_ok=True)
    d, inp = condition_inputs(geom, gcfg, obj, seed)
    n = geom.num_latents
    batch = {f"xyz_fps_{n}": inp["xyz"], f"sdf_fps_{n}": inp["sdf"], f"xyz_fps_{4 * n}": inp["xyz_kv"],
             f"sdf_fps_{4 * n}": inp["sdf_kv"]}
    code = geom_code(geom, batch)
    cond = {"code": code, "xyz": inp["xyz"], "sdf": inp["sdf"]}
    # The same condition as inference makes: domain_condition draws these sets with this seed.
    same, _, bnd = domain_condition(geom, gcfg, obj, "cpu", seed=seed)
    if not torch.equal(same["code"], code):
        raise RuntimeError("the reference's condition differs from domain_condition's: keep the two in step")
    pe = sit.cond_embed(cond)
    white = torch.randn((n, sit.channels), generator=torch.Generator().manual_seed(int(seed)))
    a, b = (torch.as_tensor(x, dtype=torch.float32) for x in bnd)
    z0 = taubin_smooth(white, inp["xyz"][0].float(), a, b)
    v0 = sit(z0[None], torch.zeros(1), cond, pe=pe)

    def f(x, t):
        return sit(x, torch.full((1,), float(t)), cond, pe=pe)

    t0 = time.time()
    x1, nfe = dopri5(f, z0[None], atol=atol, rtol=rtol)
    seconds = time.time() - t0
    latent = sit.denormalize(x1)
    again, nfe_again = sample_latents(sit, cond, [seed], "taubin", [bnd], "dopri5", atol=atol, rtol=rtol)
    if not torch.equal(again, latent) or nfe_again != nfe:
        raise RuntimeError("the reference's sample differs from sample_latents': keep the two in step")
    tokens = vae.return_last_features(latent)
    V, T, levels = extraction_mesh(obj, h_rel)
    P = ((V[T].mean(1) - d.centre) * d.scale).astype(np.float32)
    fields = decode_fields(vae, latent, P)
    log(f"  reference {to_repo(obj)}, seed {seed}: {nfe} evaluations ({seconds:.0f} s on the CPU), extraction mesh "
        f"{len(T)} triangles ({levels} subdivisions)")

    arrays = {"xyz": inp["xyz"][0], "sdf": inp["sdf"][0], "xyz_kv": inp["xyz_kv"][0], "sdf_kv": inp["sdf_kv"][0],
              "code": code[0], "pe_head": torch.stack(pe)[:, 0, :64], "noise": white, "z0": z0, "v0": v0[0],
              "x1": x1[0], "latent": latent[0], "tokens": tokens[0]}
    arrays = {k: v.float().numpy() for k, v in arrays.items()}
    arrays.update(V=np.asarray(V, np.float64), T=np.asarray(T, np.int32),
                  **{k: np.asarray(fields[k], np.float32) for k in QUERY_KEYS})
    for k, v in arrays.items():
        np.save(os.path.join(out, f"{k}.npy"), np.ascontiguousarray(v))
    shutil.copyfile(obj, os.path.join(out, "domain.obj"))

    # QUARRY on the decoded fields, as generate runs it: a subprocess, since torch is loaded here.
    work = os.path.join(out, "quarry")
    shutil.rmtree(work, ignore_errors=True)
    write_quarry_input(work, V, T, fields, d.scale, DEFAULT_SWEEPS)
    rc = run_quarry(None, work)["rc"]
    q = read_output(work)
    shutil.rmtree(work, ignore_errors=True)
    log(f"  QUARRY on the reference fields: " + (f"{q['blocks']} blocks, " + ("covers" if q["covers"] else
                                                                               "does not cover")
                                                 + f", {q['inverted']} quads inverted" if rc == 0 else f"rc {rc}"))
    case = {"obj": "domain.obj", "source": to_repo(obj), "seed": int(seed), "tokens": n,
            "centre": [float(c) for c in d.centre], "scale": float(d.scale), "h_rel": h_rel, "levels": int(levels),
            "triangles": int(len(T)), "nfe": int(nfe), "atol": atol, "rtol": rtol, "taubin": list(TAUBIN),
            "quarry_rc": int(rc), "quarry_covers": int(q["covers"]), "quarry_blocks": int(q["blocks"]),
            "quarry_failed": int(q["failed"]), "quarry_inverted": int(q["inverted"])}
    with open(os.path.join(out, "case.txt"), "w") as fh:
        fh.write("# SQuadGen reference case: python -m squadgen2d export-onnx wrote it, TestSQUADGEN --check reads it."
                 "\n# Arrays are the .npy files beside it; points and fields in the domain's normalised frame.\n")
        for k, v in case.items():
            fh.write(f"{k} " + (" ".join(repr(x) for x in v) if isinstance(v, list) else repr(v)
                                if isinstance(v, float) else str(v)) + "\n")
    return case, arrays


def check_graphs(out, case, arrays, log):
    """Each graph on the reference inputs with onnxruntime, against torch's outputs: per graph, the largest error
    over the largest value. Skipped, with a note, when onnxruntime is not installed."""
    try:
        import onnxruntime as ort
    except ImportError:
        log("onnxruntime is not installed here: the graphs are checked by TestSQUADGEN --check instead")
        return None
    sess = {g: ort.InferenceSession(os.path.join(out, f"{g}.onnx"), providers=["CPUExecutionProvider"])
            for g in GRAPHS}
    a = {k: v[None] for k, v in arrays.items() if k not in ("V", "T", "pe_head")}
    worst = {}

    def compare(name, got, want):
        err = float(np.abs(got - want).max() / max(1e-12, float(np.abs(want).max())))
        worst[name] = max(worst.get(name, 0.0), err)

    code = sess["geom_code"].run(["code"], {k: a[k] for k in ("xyz", "sdf", "xyz_kv", "sdf_kv")})[0]
    compare("geom_code", code, a["code"])
    pe = sess["sit_condition"].run(["pe"], {"code": a["code"], "xyz": a["xyz"], "sdf": a["sdf"]})[0]
    compare("sit_condition", pe[:, 0, :64], arrays["pe_head"])
    v = sess["sit_velocity"].run(["v"], {"x": a["z0"], "t": np.zeros(1, np.float32), "pe": pe})[0]
    compare("sit_velocity", v, a["v0"])
    tokens = sess["vae_tokens"].run(["tokens"], {"z": a["latent"]})[0]
    compare("vae_tokens", tokens, a["tokens"])
    V, T = arrays["V"], arrays["T"]
    P = ((V[T].mean(1) - np.asarray(case["centre"])) * case["scale"]).astype(np.float32)[None]
    for k, got in zip(QUERY_KEYS, sess["vae_query"].run(list(QUERY_KEYS), {"tokens": a["tokens"], "xy": P})):
        compare("vae_query", got, a[k])
    log(f"onnxruntime {ort.__version__} against torch (largest error / largest value): "
        + ", ".join(f"{k} {v:.1e}" for k, v in worst.items()))
    return worst


# =============================================================================
# 3. MAIN
# =============================================================================


def main(argv=None):
    ap = argparse.ArgumentParser(prog="squadgen2d export-onnx", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default=DEFAULT_MODEL, help="the SQ-Diffuse checkpoint (default best_models/best.pt)")
    ap.add_argument("--vae", help="its SQ-VAE, when not at the path the checkpoint records")
    ap.add_argument("--out", default=DEFAULT_OUT, help="where the graphs go (default best_models/onnx)")
    ap.add_argument("--opset", type=int, default=OPSET)
    ap.add_argument("--reference", default=DEFAULT_REFERENCE, help="the reference case's domain OBJ ('' for none)")
    ap.add_argument("--seed", type=int, default=0, help="the reference case's seed (token sets and noise)")
    ap.add_argument("--h_rel", type=float, default=0.002, help="the reference case's extraction mesh")
    ap.add_argument("--atol", type=float, default=1e-6)
    ap.add_argument("--rtol", type=float, default=1e-3)
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)

    def log(msg):
        print(msg, flush=True)

    t0 = time.time()
    sit, vae, geom, cfg, vcfg, gcfg = load_sq_diffuse(args.model, "cpu", args.vae)
    ck = torch.load(args.model, map_location="cpu", weights_only=False)
    vae_fp = cfg.get("vae_fingerprint") or fingerprint(vae, geom)
    meta = metadata(sit, vae, cfg, gcfg, args.model, ck["step"], vae_fp)
    log(f"export of {to_repo(args.model)} (step {ck['step']}, N = {meta['squadgen.tokens']}, SQ-VAE {cfg['vae']}) "
        f"-> {to_repo(args.out)}, opset {args.opset}")
    times = export_graphs(sit, vae, geom, meta, args.out, args.opset, log)
    summary = {"metadata": meta, "graphs": {g: {"file": f"{g}.onnx", "bytes": os.path.getsize(
        os.path.join(args.out, f"{g}.onnx")), "seconds": round(times[g], 1)} for g in GRAPHS},
        "metrics": ck.get("metrics")}
    if args.reference:
        ref = os.path.join(args.out, "reference")
        case, arrays = reference_case(sit, vae, geom, gcfg, os.path.abspath(args.reference), args.seed, args.h_rel,
                                      args.atol, args.rtol, ref, log)
        summary["reference"] = case
        worst = check_graphs(args.out, case, arrays, log)
        if worst is not None:
            summary["onnxruntime_check"] = worst
    with open(os.path.join(args.out, "squadgen.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    log(f"done in {time.time() - t0:.0f} s: {', '.join(f'{g}.onnx' for g in GRAPHS)}"
        + (", reference/" if args.reference else "") + ", squadgen.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
