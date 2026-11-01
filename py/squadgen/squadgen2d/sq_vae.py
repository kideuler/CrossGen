"""
sq_vae.py -- train SQuadGen's SQ-VAE on CrossGen layouts of planar domains.

This is the SQ-VAE of docs/squadgen_2d.md Sec. 3.3. It autoencodes a quad layout. The encoder reads the layout's chart
distance fields at N token points, conditioned on the frozen Geom-AE's geometry code at the same points, and outputs
an N x latent_dim code. The decoder rebuilds the fields at any point of Omega from that code alone. SQ-Diffuse later
learns to generate the code from the geometry, and QuadExtraction turns decoded fields back into a layout (Sec. 1.2).

The network is the release's KLAutoEncoder (SQuadGen squadgen/network/models_ae.py) as sqvae_geomae.yaml's
latent_vae_params configures it, with the Sec. 3.1 edits for 2D. Module names and shapes are the release's. The
checkpoint holds the Geom-AE under `cond_vae.` and SQ-VAE under `latent_vae.`: the state dict of the release's
KLAutoEncoderJoint, which main_ae.py saves and SQ-Diffuse loads. What is kept, and why:

* Encoder. The queries are the geometry code itself. Token i starts as proj_geom_code_fps(code_i) +
  point_embed(x_i, fields at x_i), with use_geom_code_as_fps = use_geom_code_and_pe_as_fps = 1. x_i = xyz_fps_N[i] is
  the point Geom-AE token i was computed at (Sec. 3.2, "Token alignment"). The tokens cross-attend once to the 4N
  key/value points and their fields, through one head as wide as the model, then a FeedForward, then a diagonal
  Gaussian of latent_dim channels. The code is the Geom-AE's posterior mean (use_geom_cond_mean = 1), computed
  without gradient.
* Point features (in_keys_list): dcdf, cdf, their unit gradients gdcdf and gcdf (renormalised on input,
  gcolor_grad_norm_in = 1), the seven offsets, and the truncated SDF in place of the release's normal (Sec. 3.1 item
  7). There is no context LayerNorm anywhere (is_context_norm = 0).
* Decoder: proj, `depth` pre-norm self-attention blocks (8 heads x 64), one single-head cross-attention from the query
  point's Fourier embedding, and a zero-initialised 3-layer MLP per output key (dcdf, cdf and the seven offsets). The
  decoder never sees the geometry code, so the latent has to carry the layout by itself. That is what lets
  SQ-Diffuse's samples be decoded. The offset heads divide by 7 (offset) and 3.5 (offset1, offset2), as the release's
  do (Sec. 3.3, "Offset scaling"); the release leaves offset3 and offsetb-d unscaled too.
* Loss: engine_ae.calc_loss with the Geom-AE frozen. Mean L1 per output key (per sample over points and channels,
  then over the batch), weights 1, plus kl_weight x KL with kl_weight = 1e-3. The decoder sees a posterior sample.

What changes for 2D:

* Every point is an area sample of Omega. The NPZ's fps token sets and its 20 000 sample points lie on the split quad
  mesh. Positions and every vector field are 2-d (Sec. 3.1 items 1-2).
* The frame is the Geom-AE's (area centroid at the origin, farthest boundary vertex at 0.999), not the release
  loader's vertex mean, so the frozen encoder sees domains in the frame it was trained in.
* Augmentation is the Geom-AE's (geometry.Deformation: rotation about z, reflection, anisotropic scale and a smooth
  warp, folds rejected), with every field carried along by it (cdf_pack.make_sample), as the release recomputes them
  on its add_noise output. sdf is measured to the deformed split mesh's own boundary (Appendix D).

Data: a field pack (cdf_pack.py) of a cdf-data directory: every layout's fields at its token and query points,
presaved, so that a training sample is a few array slices. A pack that is missing or no longer matches the manifest
is built first, with the Geom-AE's N and tau (--pack_jobs processes), as the Geom-AE builds its token-set cache, so
the NPZs and manifest.csv are all a machine needs (vast_sq_vae.sh uploads nothing else). Different methods give
different layouts of one domain,
and since SQ-VAE reconstructs a given layout, each is an example of its own. The split is the frozen Geom-AE's: a
layout is validation when its domain's family is one the Geom-AE held out (domains.geom_ae_val_families), so no
network in the chain trained on a validation domain. The SQ-VAE checkpoint records the families, and SQ-Diffuse
validates on the same ones. When the Geom-AE's split cannot be found, the manifest's `split` column is used, which
cdf_data.py drew the way a Geom-AE trained on the same OBJ directories would.

Run
---
    python -m squadgen2d sq-vae                                               # train -> py/squadgen/runs/sq_vae_2d
    python -m squadgen2d cdf-pack py/squadgen/data/cdf                        # the pack on its own (sq-vae builds it)
    python -m squadgen2d sq-vae --data py/squadgen/data/cdf_oracle --geom_ae py/squadgen/runs/geom_ae_2d_current/best.pt \\
        --out py/squadgen/runs/sq_vae_2d_oracle2 --steps 30000
    python -m squadgen2d sq-vae --out /tmp/smoke --steps 200                  # smoke test
    python -m squadgen2d sq-vae --eval_only                                   # metrics, plots, latent stats of best.pt
    python -m squadgen2d sq-vae --data py/squadgen/data/cdf_oracle --init py/squadgen/runs/sq_vae_2d_oracle/best.pt \\
        --out py/squadgen/runs/sq_vae_2d_oracle_more --steps 10000 --lr 1e-4    # more steps from a checkpoint

Continuing (--init), as the Geom-AE does. A checkpoint holds the weights but no optimizer state, so a continuation is
a new run that starts from them: a fresh AdamW with its own warm-up and cosine over --steps (step numbers start
again at 1; the config records init and init_step). The network's shape comes from the checkpoint, and so does the
frozen Geom-AE: the one inside it, which these weights were trained against, whatever --geom_ae says. So does the
split, when the checkpoint records its held-out families. --data may differ (more layouts, say): the split is by the
Geom-AE's families, so it cannot leak a validation domain into training. The loaded model is evaluated at step 0
and saved as the first best.pt, so best.pt never ends worse than the start. This is also how to carry on a cloud run
that was cut off: its last best.pt or last.pt is the start.

Outputs in --out:

* best.pt / last.pt: checkpoint with the joint state dict, both configs, step and metrics;
* recon_val.png / recon_train.png: true and decoded CDF and DCDF per layout;
* latent_stats.json: SQ-Diffuse's target_mu / target_std (Sec. 3.4);
* log.txt.

Memory and time on the 8 GB M2 (MPS), batch 8: ~0.8 s/step at n_query 2048 (3.3 GiB), 20k steps ~4.3 h. Downstream
code calls load_sq_vae(), encode_layout() and decode_fields().
"""

import argparse
import json
import math
import os
import time

import numpy as np
import torch
import torch.nn as nn

from . import training as tr
from .cdf_pack import QUERY_FIELDS, default_pack, make_sample, open_pack, record_from_npz
from .domains import default_fps_workers, geom_ae_val_families
from .geom_ae import build_model as build_geom
from .geom_ae import load_geom_ae
from .layers import KEY_DIM_DICT, Attention, FeedForward, PointEmbed, PreNorm, xavier_init
from .layout import METHODS, OFFSET_KEYS, Layout, chart_sides, layout_fields
from .paths import DEFAULT_CDF, DEFAULT_GEOM_AE, DEFAULT_SQ_VAE, from_repo, portable, to_repo

IN_KEYS = ("dcdf", "gdcdf", "cdf", "gcdf", "sdf") + OFFSET_KEYS
OUT_KEYS = ("dcdf", "cdf", "offset", "offset1", "offset2", "offset3", "offsetb", "offsetc", "offsetd")
assert set(OUT_KEYS) == {k for k, _ in QUERY_FIELDS} - {"xyz"}, "the pack's query fields are SQ-VAE's outputs"
EVAL_QUERIES = 4096


# =============================================================================
# 1. DATA
# =============================================================================


class LayoutStream(torch.utils.data.IterableDataset):
    """Endless stream of samples drawn uniformly over the layouts `indices` of a field pack, augmented when
    cfg["augment"]. Each DataLoader worker seeds its own generator from (seed, worker id)."""

    def __init__(self, pack, indices, cfg):
        self.pack, self.indices, self.cfg = pack, list(indices), cfg

    def __iter__(self):
        info = torch.utils.data.get_worker_info()
        rng = np.random.default_rng([self.cfg["seed"], info.id if info is not None else 0])
        while True:
            i = self.indices[rng.integers(len(self.indices))]
            yield make_sample(self.pack.record(i), self.cfg, rng, self.cfg["augment"])


def split_layouts(rows, val_fams):
    """(train, val) row indices: val where the domain's family is held out; the manifest's `split` without one."""
    is_val = (lambda r: r["family"] in val_fams) if val_fams is not None else (lambda r: r["split"] == "val")
    return [i for i, r in enumerate(rows) if not is_val(r)], [i for i, r in enumerate(rows) if is_val(r)]


# =============================================================================
# 2. NETWORK  (the release's KLAutoEncoder as SQ-VAE, models_ae.py, with the Sec. 3.1 edits)
# =============================================================================


class SQVAE(nn.Module):
    """The release's KLAutoEncoder with sqvae_geomae.yaml's latent_vae_params: is_context_norm = 0 (so no PreNorm
    normalises a context; the decoder still normalises its queries), use_geom_code_as_fps =
    use_geom_code_and_pe_as_fps = 1, gcolor_grad_norm_in = 1. No learnable latents, no encoder self-attention, no
    decoder FeedForward, no QK-norm. proj_geom_code_fps maps geom_dim to dim; the release has both 512, and here both
    are 256."""

    def __init__(self, in_keys_list=IN_KEYS, out_keys_list=OUT_KEYS, out_weights_list=None, depth=8, dim=256,
                 geom_dim=256, num_latents=512, latent_dim=16, heads=8, dim_head=64, drop_path_rate=0.1):
        super().__init__()
        self.num_latents = num_latents
        self.in_keys_list_mlp = sorted(set(in_keys_list))
        self.out_keys_list = list(out_keys_list)
        self.out_weights_list = list(out_weights_list) if out_weights_list else [1.0] * len(self.out_keys_list)
        self.out_keys_list_mlp = sorted(set(out_keys_list))

        self.cross_attend_blocks = nn.ModuleList([
            PreNorm(dim, Attention(dim, dim, heads=1, dim_head=dim), context_dim=dim, is_context_norm=0),
            PreNorm(dim, FeedForward(dim)),
        ])
        self.point_embed = PointEmbed(dim, self.in_keys_list_mlp)
        self.layers = nn.ModuleList([nn.ModuleList([
            PreNorm(dim, Attention(dim, heads=heads, dim_head=dim_head, drop_path_rate=drop_path_rate)),
            PreNorm(dim, FeedForward(dim, drop_path_rate=drop_path_rate)),
        ]) for _ in range(depth)])
        self.decoder_cross_attn = PreNorm(dim, Attention(dim, dim, heads=1, dim_head=dim), context_dim=dim,
                                          is_context_norm=0)
        self.to_outputs = nn.ModuleDict({k: nn.Sequential(
            nn.Linear(dim, dim), nn.ReLU(), nn.Linear(dim, dim), nn.ReLU(), nn.Linear(dim, KEY_DIM_DICT[k]))
            for k in self.out_keys_list_mlp})
        self.proj = nn.Linear(latent_dim, dim)
        self.mean_fc = nn.Linear(dim, latent_dim)
        self.logvar_fc = nn.Linear(dim, latent_dim)
        self.proj_geom_code_fps = nn.Linear(geom_dim, dim)
        self.initialize_weights()

    def initialize_weights(self):
        xavier_init(self)
        nn.init.normal_(self.point_embed.mlp_pos.weight, std=1)
        for k in self.in_keys_list_mlp:
            nn.init.normal_(getattr(self.point_embed, f"mlp_{k}")[0].weight, std=0.02)
        for k in self.out_keys_list_mlp:
            nn.init.constant_(self.to_outputs[k][-1].weight, 0)
            nn.init.constant_(self.to_outputs[k][-1].bias, 0)

    def features(self, batch, n):
        """The point features of one fps set, gradients renormalised (gcolor_grad_norm_in)."""
        out = {}
        for k in self.in_keys_list_mlp:
            x = batch[f"{k}_fps_{n}"]
            if k in ("gcdf", "gdcdf"):
                x = x / (x.norm(dim=-1, keepdim=True) + 1e-6)
            out[k] = x
        return out

    def encode(self, batch, geom_code, is_mode=False):
        """-> {"latent": (B, N, latent_dim), "kl": (B,), "mean", "logvar"}. geom_code is the frozen Geom-AE's
        (B, N, geom_dim) code at xyz_fps_N (geom_code())."""
        n = self.num_latents
        x = self.proj_geom_code_fps(geom_code) + self.point_embed(batch[f"xyz_fps_{n}"], self.features(batch, n))
        pc = self.point_embed(batch[f"xyz_fps_{4 * n}"], self.features(batch, 4 * n))
        cross_attn, cross_ff = self.cross_attend_blocks
        x = cross_attn(x, context=pc) + x
        x = cross_ff(x) + x
        mean, logvar = self.mean_fc(x), self.logvar_fc(x).clamp(-10.0, 10.0)
        z = mean if is_mode else mean + (0.5 * logvar).exp() * torch.randn_like(mean)
        kl = 0.5 * (mean ** 2 + logvar.exp() - 1.0 - logvar).mean(dim=[1, 2])
        return {"latent": z, "kl": torch.nan_to_num(kl, 0.0, 0.0, 0.0), "mean": mean, "logvar": logvar}

    def return_last_features(self, z):
        x = self.proj(z)
        for self_attn, self_ff in self.layers:
            x = self_attn(x) + x
            x = self_ff(x) + x
        return x

    def query_last_features(self, x, queries):
        h = self.decoder_cross_attn(self.point_embed.query(queries), context=x)
        out = {}
        for k in self.out_keys_list_mlp:
            y = self.to_outputs[k](h)
            if k == "offset":
                y = y / 7
            elif k in ("offset1", "offset2"):
                y = y / 3.5
            out[k] = y
        return out

    def decode(self, z, queries):
        """Every output key at queries (B, P, 2): what QuadExtraction will read off a triangulation of Omega (Sec. 4
        step 6)."""
        return self.query_last_features(self.return_last_features(z), queries)


def geom_code(geom, batch):
    """The frozen Geom-AE's code at xyz_fps_N from its posterior mean (use_geom_cond_mean=1)."""
    with torch.no_grad():
        return geom.return_last_features_and_norm(geom.encode(batch, is_mode=True)["latent"])


MODEL_KEYS = ("depth", "dim", "geom_dim", "num_latents", "latent_dim", "heads", "dim_head", "drop_path_rate")
# The MODEL_KEYS that are SQ-VAE's own options (geom_dim and num_latents are the Geom-AE's): what --init takes from
# the checkpoint, since its weights fit only that network.
ARCH_KEYS = ("depth", "dim", "latent_dim", "heads", "dim_head", "drop_path_rate")


def build_model(cfg):
    return SQVAE(**{k: cfg[k] for k in MODEL_KEYS})


def save(path, model, geom, cfg, gcfg, step, metrics):
    """main_ae.py's layout: KLAutoEncoderJoint's state dict, the Geom-AE under `cond_vae.` and SQ-VAE under
    `latent_vae.`."""
    sd = {f"cond_vae.{k}": v.detach().cpu() for k, v in geom.state_dict().items()}
    sd.update({f"latent_vae.{k}": v.detach().cpu() for k, v in model.state_dict().items()})
    torch.save({"model": sd, "config": portable(cfg), "geom_config": portable(gcfg), "step": step,
                "metrics": metrics}, path)


def load_sq_vae(path, device="cpu"):
    """SQ-VAE and its frozen Geom-AE from one checkpoint (its path, or the dict torch.load returned): returns (model,
    geom, config)."""
    ck = path if isinstance(path, dict) else torch.load(path, map_location=device, weights_only=False)
    geom = build_geom(ck["geom_config"]).to(device)
    geom.load_state_dict({k[len("cond_vae."):]: v for k, v in ck["model"].items() if k.startswith("cond_vae.")})
    geom.requires_grad_(False)
    model = build_model(ck["config"]).to(device)
    model.load_state_dict({k[len("latent_vae."):]: v for k, v in ck["model"].items() if k.startswith("latent_vae.")})
    return model.eval(), geom.eval(), ck["config"]


@torch.no_grad()
def encode_layout(model, geom, cfg, npz_path, device="cpu", seed=0):
    """The posterior-mean latent (1, N, latent_dim) of one layout NPZ, unaugmented, in the Geom-AE's frame. Returns
    (latent, batch, Layout); Layout.centre/scale map the frame back."""
    rec, L = record_from_npz(npz_path, cfg["n_tokens"], cfg["tau"])
    s = make_sample(rec, dict(cfg, n_query=1), np.random.default_rng(seed), False)
    b = tr.to_dev({k: v[None] for k, v in s.items()}, device)
    return model.encode(b, geom_code(geom, b), is_mode=True)["latent"], b, L


@torch.no_grad()
def decode_fields(model, z, xy, chunk=8192):
    """Decode a latent (1, N, latent_dim) at points xy (P, 2) of the normalised frame, in chunks: {key: (P, channels)
    numpy}."""
    dev = z.device
    parts = [model.decode(z, torch.as_tensor(xy[s:s + chunk], dtype=torch.float32, device=dev)[None])
             for s in range(0, len(xy), chunk)]
    return {k: torch.cat([p[k][0] for p in parts]).cpu().numpy() for k in parts[0]}


# =============================================================================
# 3. EVALUATION
# =============================================================================


def calc_loss(model, geom, b, kl_weight, is_mode=False):
    """engine_ae.calc_loss with a frozen Geom-AE: mean L1 per output key (per sample over points and channels, then
    over the batch), weighted, plus kl_weight x KL."""
    enc = model.encode(b, geom_code(geom, b), is_mode=is_mode)
    pred = model.decode(enc["latent"], b["xyz_query"])
    terms = {k: (pred[k] - b[f"{k}_query"]).abs().mean(dim=[1, 2]).mean() for k in model.out_keys_list}
    terms["kl"] = enc["kl"].mean()
    loss = sum(w * terms[k] for k, w in zip(model.out_keys_list, model.out_weights_list))
    return loss + kl_weight * terms["kl"], terms


class LayoutEval:
    """Fixed, unaugmented samples of some of a pack's layouts, decoded from the posterior mean at EVAL_QUERIES area
    samples each. Per layout:
      l1        weighted mean over the output keys of their L1 (best.pt is chosen on its mean)
      l1_<key>  each key's L1: dcdf and cdf are in [0, 1], offsets in the frame's units
      wall_iou  IoU of the decoded walls {cdf <= wall} with the true ones. wall = 0.1 is where QuadExtraction's region
                growing stops (-t 0.1, Sec. 5.1 step 2), so this is how well the decoded chart sides sit where the
                true ones are."""

    def __init__(self, pack, indices, cfg):
        self.rows, self.items, self.wall = [pack.rows[i] for i in indices], [], cfg["wall"]
        small = dict(cfg, n_query=EVAL_QUERIES)
        for k, i in enumerate(indices):
            self.items.append(make_sample(pack.record(i), small, np.random.default_rng(1000 + k), False))

    @torch.no_grad()
    def __call__(self, model, geom, dev, keep=False):
        was = model.training
        model.eval()
        per, latents = [], []
        w = dict(zip(model.out_keys_list, model.out_weights_list))
        for s in self.items:
            b = tr.to_dev({k: v[None] for k, v in s.items()}, dev)
            z = model.encode(b, geom_code(geom, b), is_mode=True)["latent"]
            pred = model.decode(z, b["xyz_query"])
            m = {f"l1_{k}": (pred[k] - b[f"{k}_query"]).abs().mean().item() for k in model.out_keys_list}
            m["l1"] = sum(w[k] * m[f"l1_{k}"] for k in w) / sum(w.values())
            pw, tw = pred["cdf"][0, :, 0].cpu().numpy() <= self.wall, s["cdf_query"][:, 0] <= self.wall
            m["wall_iou"] = float((pw & tw).sum() / max(1, (pw | tw).sum()))
            per.append(m)
            if keep:
                latents.append(z.cpu())
        model.train(was)
        res = {k: float(np.mean([m[k] for m in per])) for k in per[0]}
        return (res, per, latents) if keep else res


# Sequential blue, field value 0 (light) to 1 (dark): the dataviz skill's default ramp.
FIELD_RAMP = ("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")


def plot_reconstructions(ev, per, latents, model, dev, out, title, k=6):
    """A row per layout, worst L1 first, so that the failures are the first thing seen; one layout per domain, since
    the methods often agree on a domain's layout. The columns are the true and decoded CDF and DCDF, coloured per
    triangle of the split mesh at its centroid, with the true chart sides drawn over each. A row is as tall as its
    domain's aspect. Needs the layouts' NPZs (the split meshes), not just the pack."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap

    order = list(np.argsort([m["l1"] for m in per])[::-1])
    seen, worst, best = set(), [], []
    for pool, n, picked in ((order, k // 2, worst), (order[::-1], k - k // 2, best)):
        for i in pool:
            if len(picked) < n and ev.rows[i]["model"] not in seen:
                seen.add(ev.rows[i]["model"])
                picked.append(i)
    idx = worst + best[::-1]
    layouts = [Layout(ev.rows[i]["npz"]) for i in idx]
    aspect = [float(np.clip(np.ptp(L.V[:, 1]) / np.ptp(L.V[:, 0]), 0.3, 1.3)) for L in layouts]
    cmap = LinearSegmentedColormap.from_list("field", FIELD_RAMP)
    cols = ("true CDF", "decoded CDF", "true DCDF", "decoded DCDF")
    fig, axes = plt.subplots(len(idx), 4, figsize=(10.4, sum(2.4 * a + 0.45 for a in aspect) + 0.6), squeeze=False,
                             gridspec_kw={"height_ratios": [a + 0.19 for a in aspect]})
    im = None
    for row, i, L in zip(axes, idx, layouts):
        r, m = ev.rows[i], per[i]
        tri = np.arange(len(L.T))
        true = layout_fields(L, L.V, tri, np.full((len(tri), 3), 1 / 3))
        dec = decode_fields(model, latents[i].to(dev), true["xyz"])
        sides = chart_sides(r["npz"], L.V)
        for ax, vals in zip(row, (true["cdf"], dec["cdf"], true["dcdf"], dec["dcdf"])):
            im = ax.tripcolor(L.V[:, 0], L.V[:, 1], L.T, facecolors=np.clip(vals[:, 0], 0, 1), cmap=cmap,
                              vmin=0, vmax=1, shading="flat", edgecolors="none", rasterized=True)
            ax.add_collection(LineCollection(sides, colors="#1a1a19", linewidths=0.6))
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color("#c9c8c4")
        charts = f"{r['charts']} chart" + ("" if r["charts"] == "1" else "s")
        row[0].set_ylabel(f"{r['model']} {r['method']}\n{charts}", fontsize=8, color="#1a1a19")
        row[1].set_title(f"L1 {m['l1']:.3f}   wall IoU {m['wall_iou']:.3f}", fontsize=8, color="#52514e")
    for ax, c in zip(axes[0], cols):
        ax.set_title(c + ("\n" + ax.get_title() if ax.get_title() else ""), fontsize=8, color="#1a1a19")
    fig.suptitle(f"{title}\nlines: true chart sides", y=0.997, fontsize=9, color="#1a1a19")
    fig.tight_layout(rect=(0, 0, 0.92, 0.975))
    cax = fig.add_axes((0.935, 0.25, 0.012, 0.5))
    cb = fig.colorbar(im, cax=cax, ticks=[0, 0.25, 0.5, 0.75, 1])
    cb.set_label("field value", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.savefig(out, dpi=110, facecolor="#fcfcfb")
    plt.close(fig)


@torch.no_grad()
def latent_stats(model, geom, pack, indices, cfg, dev, n):
    """SQ-Diffuse's target_mu / target_std (Sec. 3.4): the per-channel mean and std of the posterior mean over n
    samples of the training layouts, augmented as SQ-Diffuse will see them. The release's values are for its 3D
    SQ-VAE and must not be reused."""
    model.eval()
    rng = np.random.default_rng(cfg["seed"] + 1)
    s1 = s2 = 0.0
    count = 0
    for _ in range(math.ceil(n / cfg["batch"])):
        samples = [make_sample(pack.record(indices[rng.integers(len(indices))]), dict(cfg, n_query=1), rng,
                               cfg["augment"]) for _ in range(cfg["batch"])]
        b = tr.to_dev({k: np.stack([s[k] for s in samples]) for k in samples[0]}, dev)
        z = model.encode(b, geom_code(geom, b), is_mode=True)["latent"].cpu().double().flatten(0, 1)
        s1, s2, count = s1 + z.sum(0), s2 + (z * z).sum(0), count + len(z)
    mu = s1 / count
    std = (s2 / count - mu * mu).clamp_min(0).sqrt()
    return {"target_mu": mu.cpu().tolist(), "target_std": std.cpu().tolist(), "samples": n, "tokens": count}


# =============================================================================
# 4. TRAINING
# =============================================================================


def eval_line(step, mv, mt):
    return (f"  eval {step:6d}  val l1 {mv['l1']:.4f} cdf {mv['l1_cdf']:.4f} dcdf {mv['l1_dcdf']:.4f} "
            f"wall_iou {mv['wall_iou']:.4f}  |  train l1 {mt['l1']:.4f} cdf {mt['l1_cdf']:.4f} "
            f"dcdf {mt['l1_dcdf']:.4f} wall_iou {mt['wall_iou']:.4f}")


def train(cfg, gcfg, geom, pack, train_idx, val_eval, train_eval, log, init=None):
    """init: an SQ-VAE checkpoint dict to continue from (see Continuing in the module docstring); geom is then the
    Geom-AE inside it."""
    dev = cfg["device"]
    amp = tr.amp_dtype(dev, cfg["amp"])
    torch.manual_seed(cfg["seed"])
    model = build_model(cfg).to(dev)
    geom = geom.to(dev)
    log(f"device={dev}  amp={amp}  workers={cfg['workers']}  SQ-VAE parameters="
        f"{sum(p.numel() for p in model.parameters()):,}  (Geom-AE frozen: {sum(p.numel() for p in geom.parameters()):,})")
    best = float("inf")
    if init is not None:
        model.load_state_dict({k[len("latent_vae."):]: v for k, v in init["model"].items()
                               if k.startswith("latent_vae.")})
        # The start is the best so far, so a continuation that only gets worse leaves best.pt as it began.
        mv, mt = val_eval(model, geom, dev), train_eval(model, geom, dev)
        log(eval_line(0, mv, mt) + f"  (start: {cfg['init']}, step {init['step']} of its run)")
        best = mv["l1"]
        save(os.path.join(cfg["out"], "best.pt"), model, geom, cfg, gcfg, 0, {"val": mv, "train": mt})
        tr.empty_cache(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    scaler = tr.grad_scaler(dev, amp)
    loader = tr.make_loader(LayoutStream(pack, train_idx, cfg), cfg)
    t0, run = time.time(), {}
    for step, batch in enumerate(loader, 1):
        for g in opt.param_groups:
            g["lr"] = tr.lr_at(step - 1, cfg)
        # The frozen Geom-AE runs under the same autocast, as the release runs both halves.
        with tr.autocast(dev, amp):
            loss, terms = calc_loss(model, geom, tr.to_dev(batch, dev), cfg["kl_weight"])
        gnorm = tr.train_step(model, opt, scaler, loss, cfg["clip_grad"])
        # The log line keeps the two fields apart and folds the seven offsets into one mean.
        shown = {"loss": loss, "dcdf": terms["dcdf"], "cdf": terms["cdf"],
                 "offsets": sum(terms[k] for k in OFFSET_KEYS) / len(OFFSET_KEYS), "kl": terms["kl"], "gnorm": gnorm}
        tr.accumulate(run, shown)
        if step % cfg["log_every"] == 0 or step == 1:
            el = time.time() - t0
            log(f"step {step:6d}  lr {opt.param_groups[0]['lr']:.2e}  {tr.averages(run)}"
                + f"  ({el / step * 1000:.0f} ms/step, {el / 60:.1f} min{tr.memory_note(dev)})")
            run = {}
        if step % cfg["eval_every"] == 0 or step == cfg["steps"]:
            mv, mt = val_eval(model, geom, dev), train_eval(model, geom, dev)
            log(eval_line(step, mv, mt))
            save(os.path.join(cfg["out"], "last.pt"), model, geom, cfg, gcfg, step, {"val": mv, "train": mt})
            if mv["l1"] < best:
                best = mv["l1"]
                save(os.path.join(cfg["out"], "best.pt"), model, geom, cfg, gcfg, step, {"val": mv, "train": mt})
            tr.empty_cache(dev)
        if step >= cfg["steps"]:
            break


def main(argv=None):
    ap = argparse.ArgumentParser(prog="squadgen2d sq-vae", description=__doc__.split("\n")[1].strip())
    ap.add_argument("--data", default=DEFAULT_CDF, help="a cdf-data output directory (manifest.csv + npz/)")
    ap.add_argument("--pack", help="its field pack (default: DATA/pack_N<the Geom-AE's N>; see cdf-pack)")
    ap.add_argument("--geom_ae", help="Geom-AE checkpoint, frozen (default: py/squadgen/runs/geom_ae_2d/best.pt; "
                                      "with --init or --eval_only, the one inside the SQ-VAE checkpoint)")
    ap.add_argument("--out", default=DEFAULT_SQ_VAE)
    ap.add_argument("--methods", nargs="+", default=list(METHODS), choices=METHODS)
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--batch", type=int, default=8)                  # release: 16 per GPU
    ap.add_argument("--n_query", type=int, default=2048)             # release: 8192
    ap.add_argument("--dim", type=int, default=256)                  # Sec. 3.5: 256 / 8 / 16
    ap.add_argument("--depth", type=int, default=8)
    ap.add_argument("--latent_dim", type=int, default=16)
    ap.add_argument("--heads", type=int, default=8)                  # release hard-codes 8 x 64
    ap.add_argument("--dim_head", type=int, default=64)
    ap.add_argument("--drop_path_rate", type=float, default=0.1)     # release default
    ap.add_argument("--kl_weight", type=float, default=1e-3)         # sqvae_geomae.yaml
    ap.add_argument("--lr", type=float, default=2e-4)                # release: 1e-4
    ap.add_argument("--min_lr", type=float, default=1e-6)
    ap.add_argument("--warmup_frac", type=float, default=0.03)
    ap.add_argument("--wd", type=float, default=0.01)                # torch AdamW default, as the release
    ap.add_argument("--clip_grad", type=float, default=1.0)
    ap.add_argument("--no_aug", action="store_true")
    ap.add_argument("--wall", type=float, default=0.1, help="CDF wall threshold of the eval's wall IoU (QuadExtraction -t)")
    ap.add_argument("--stats_samples", type=int, default=1024, help="samples for latent_stats.json")
    ap.add_argument("--pack_jobs", type=int,
                    help="processes when the pack has to be built (default: every core on CUDA, 6 elsewhere)")
    tr.add_device_args(ap)
    ap.add_argument("--eval_only", action="store_true", help="evaluate, plot and write latent stats of best.pt")
    ap.add_argument("--init", help="continue from this SQ-VAE checkpoint (best.pt / last.pt; see Continuing)")
    args = ap.parse_args(argv)

    best = os.path.join(args.out, "best.pt")
    init, notes = None, []
    if args.eval_only or args.init:
        # Evaluating a run (perhaps trained elsewhere and copied here), or continuing one: the Geom-AE inside it,
        # which its weights were trained against, and its split.
        if args.eval_only:
            src = best
        else:
            src = args.init = args.init if os.path.exists(args.init) else from_repo(args.init)
        if not os.path.exists(src):
            raise SystemExit(f"no SQ-VAE checkpoint {src}" + (": train one first, or pass --out" if args.eval_only
                                                               else ""))
        ck = torch.load(src, map_location="cpu", weights_only=False)
        gcfg, geom = ck["geom_config"], None
        recorded = ck["config"].get("val_families")
        if args.geom_ae:
            notes.append(f"note: --geom_ae {args.geom_ae} ignored: the Geom-AE is the one inside {src}")
        args.geom_ae = ck["config"].get("geom_ae")
        if args.init and not args.eval_only:
            init = ck
            _, geom, _ = load_sq_vae(ck)
            for k in ARCH_KEYS:
                setattr(args, k, ck["config"][k])   # the network the weights fit
    else:
        args.geom_ae = args.geom_ae or os.path.join(DEFAULT_GEOM_AE, "best.pt")
        args.geom_ae = args.geom_ae if os.path.exists(args.geom_ae) else from_repo(args.geom_ae)
        if not os.path.exists(args.geom_ae):
            raise SystemExit(f"no Geom-AE checkpoint {args.geom_ae}: train one with  python -m squadgen2d geom-ae")
        geom, gcfg = load_geom_ae(args.geom_ae)
        recorded = None
    dev = tr.pick_device(args.device)
    # Token i of SQ-VAE reads the point Geom-AE token i was computed at, so N is the Geom-AE's (Sec. 3.2, "Token
    # alignment"), and sdf is truncated at its tau.
    n = gcfg["num_latents"]
    cfg = dict(vars(args), pack=os.path.abspath(args.pack or default_pack(args.data, n)), n_tokens=n, num_latents=n,
               geom_dim=gcfg["dim"], tau=gcfg["tau"], augment=not args.no_aug, device=dev,
               workers=tr.default_workers(dev) if args.workers is None else args.workers,
               init_step=init["step"] if init else None)

    log = tr.RunLog(args.out)
    log.header(cfg)
    for s in notes:
        log(s)
    # Built here when it is missing or stale: a one-off CPU job before any model runs, so it may use more processes
    # than the DataLoader, as the Geom-AE's token-set cache does.
    pack = open_pack(args.data, args.pack, n, gcfg["tau"],
                     args.pack_jobs or default_fps_workers(dev), log)
    if recorded is not None:
        val_fams, how = set(recorded), "recorded in the SQ-VAE checkpoint"
    else:
        val_fams, how = geom_ae_val_families(gcfg, None if init or args.eval_only else args.geom_ae)
    keep = [i for i, r in enumerate(pack.rows) if r["method"] in args.methods]
    tr_i, va_i = split_layouts([pack.rows[i] for i in keep], val_fams)
    train_idx, val_idx = [keep[i] for i in tr_i], [keep[i] for i in va_i]
    if val_fams is None:
        log(f"warning: the Geom-AE's held-out families are unknown ({how}); using the manifest's split, which may "
            "validate on domains the Geom-AE trained on")
    else:
        log(f"split: the Geom-AE's held-out families, {how}")
    cfg["val_families"] = sorted({pack.rows[i]["family"] for i in val_idx}) if val_fams is None else sorted(val_fams)
    rows = [pack.rows[i] for i in keep]
    log(f"{len(keep)} layouts of {len({r['model'] for r in rows})} domains in {to_repo(pack.path)}: {len(train_idx)} "
        f"train, {len(val_idx)} val ({' '.join(sorted({pack.rows[i]['family'] for i in val_idx}))})")
    if not val_idx:
        raise SystemExit("no validation layouts: none of these domains is in a family the Geom-AE held out")
    rng = np.random.default_rng(args.seed)
    val_eval = LayoutEval(pack, val_idx, cfg)
    train_eval = LayoutEval(pack, [train_idx[i] for i in rng.choice(len(train_idx), min(48, len(train_idx)),
                                                                    replace=False)], cfg)

    if not args.eval_only:
        train(cfg, gcfg, geom, pack, train_idx, val_eval, train_eval, log, init)
    model, geom, _ = load_sq_vae(best, dev)
    step = torch.load(best, map_location="cpu", weights_only=False)["step"]
    # SQ-Diffuse cannot train without these, and the plots below are the likelier part to fail, so they come first.
    stats = latent_stats(model, geom, pack, train_idx, cfg, dev, args.stats_samples)
    with open(os.path.join(args.out, "latent_stats.json"), "w") as f:
        json.dump(dict(stats, step=step), f, indent=1)
    log(f"latent stats over {stats['tokens']} tokens: |mu| max {np.abs(stats['target_mu']).max():.3f}, "
        f"std {np.min(stats['target_std']):.3f}..{np.max(stats['target_std']):.3f} -> latent_stats.json")
    plots = all(os.path.exists(r["npz"]) for r in rows)
    for name, ev in (("val", val_eval), ("train", train_eval)):
        m, per, latents = ev(model, geom, dev, keep=True)
        log(f"best.pt (step {step}) on {name}: l1 {m['l1']:.4f}  wall_iou {m['wall_iou']:.4f}  "
            + "  ".join(f"{k} {m[f'l1_{k}']:.4f}" for k in OUT_KEYS)
            + f"  (wall IoU min {np.min([p['wall_iou'] for p in per]):.3f}, "
              f"median {np.median([p['wall_iou'] for p in per]):.3f})")
        if plots:
            plot_reconstructions(ev, per, latents, model, dev, os.path.join(args.out, f"recon_{name}.png"),
                                 f"SQ-VAE 2D, best.pt (step {step}), {name} layouts: worst and best L1")
    if not plots:
        log(f"no plots: the pack's NPZs are not all in {pack.data} (the plots draw their split meshes)")
    log(f"done: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
