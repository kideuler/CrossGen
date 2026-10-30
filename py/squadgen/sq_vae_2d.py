"""
sq_vae_2d.py -- train SQuadGen's SQ-VAE on CrossGen layouts of planar domains.

This is the SQ-VAE of docs/squadgen_2d.md Sec. 3.3. It autoencodes a quad layout. The
encoder reads the layout's chart distance fields at N token points, conditioned on the frozen
Geom-AE's geometry code at the same points, and outputs an N x latent_dim code. The decoder
rebuilds the fields at any point of Omega from that code alone. SQ-Diffuse later learns to
generate the code from the geometry, and QuadExtraction turns decoded fields back into a
layout (Sec. 1.2).

The training data are the NPZs that make_cdf_data.sh writes. cdf_data.py takes every
method's block decomposition of every singlemat + mambo model, with the blocks as charts, and
BlockDecomposition.cdf() builds CDFGen's split mesh and coordinate pairs in C++. The geometry
encoder is the Geom-AE that geom_ae_2d.py trained.

The network is the release's KLAutoEncoder (SQuadGen squadgen/network/models_ae.py) as
sqvae_geomae.yaml's latent_vae_params configures it, with the Sec. 3.1 edits for 2D. Module
names and shapes are the release's. The checkpoint holds the Geom-AE under `cond_vae.` and
SQ-VAE under `latent_vae.`: the state dict of the release's KLAutoEncoderJoint, which
main_ae.py saves and SQ-Diffuse loads. What is kept, and why:

* Encoder. The queries are the geometry code itself. Token i starts as
  proj_geom_code_fps(code_i) + point_embed(x_i, fields at x_i), with
  use_geom_code_as_fps = use_geom_code_and_pe_as_fps = 1. x_i = xyz_fps_N[i] is the point
  Geom-AE token i was computed at (Sec. 3.2, "Token alignment"). The tokens cross-attend once
  to the 4N key/value points and their fields, through one head as wide as the model, then a
  FeedForward, then a diagonal Gaussian of latent_dim channels. The code is the Geom-AE's
  posterior mean (use_geom_cond_mean = 1), computed without gradient.
* Point features (in_keys_list): dcdf, cdf, their unit gradients gdcdf and gcdf
  (renormalised on input, gcolor_grad_norm_in = 1), the seven offsets, and the truncated SDF
  in place of the release's normal (Sec. 3.1 item 7). There is no context LayerNorm anywhere
  (is_context_norm = 0).
* Decoder: proj, `depth` pre-norm self-attention blocks (8 heads x 64), one single-head
  cross-attention from the query point's Fourier embedding, and a zero-initialised 3-layer
  MLP per output key (dcdf, cdf and the seven offsets). The decoder never sees the geometry
  code, so the latent has to carry the layout by itself. That is what lets SQ-Diffuse's
  samples be decoded. The offset heads divide by 7 (offset) and 3.5 (offset1, offset2), as
  the release's do (Sec. 3.3, "Offset scaling"); the release leaves offset3 and offsetb-d
  unscaled too.
* Loss: engine_ae.calc_loss with the Geom-AE frozen. Mean L1 per output key (per sample
  over points and channels, then over the batch), weights 1, plus kl_weight x KL with
  kl_weight = 1e-3. The decoder sees a posterior sample.

What changes for 2D:

* Every point is an area sample of Omega. The NPZ's fps_*_{N} / fps_*_{4N} token sets and
  its 20 000 sample_* points lie on the split quad mesh. Positions and every vector field are
  2-d (Sec. 3.1 items 1-2).
* The frame is the Geom-AE's (area centroid at the origin, farthest boundary vertex at
  0.999), not the release loader's vertex mean, so the frozen encoder sees domains in the
  frame it was trained in.
* Augmentation is geom_ae_2d's (rotation about z, reflection, anisotropic scale and a smooth
  warp, with folds rejected), applied to the split mesh. Every field is then recomputed on the
  deformed mesh (Appendix B), as the release recomputes them on its add_noise output. sdf is
  measured to the deformed split mesh's own boundary (Appendix D).

Data: the status-ok rows of py/squadgen/data/cdf/manifest.csv. Train and val come from its
`split` column, which is the Geom-AE's split by family, so the Geom-AE never trained on a val
layout's domain either. Different methods give different layouts of one domain, and since
SQ-VAE reconstructs a given layout, each is an example of its own. Every sample reads its NPZ
from disk and nothing is cached, because the corpus (~2 GB) would not fit beside the model in
8 GB once per DataLoader worker.

The release trainer needs CUDA and torchrun (Sec. 3.4, "Hardware"), so this script is a plain
loop around calc_loss, as the guide suggests. It runs on CUDA, Apple MPS or the CPU (--device),
under bf16 autocast on a CUDA GPU that has it (--amp; geom_ae_2d.amp_dtype says why not on MPS),
and its checkpoints load on any of them: weights on the CPU, paths relative to the repository.

Run
---
    python py/squadgen/sq_vae_2d.py                                   # train -> py/squadgen/runs/sq_vae_2d
    python py/squadgen/sq_vae_2d.py --out /tmp/smoke --steps 200      # smoke test
    python py/squadgen/sq_vae_2d.py --eval_only                       # metrics, plots, latent stats of best.pt

Outputs in --out:

* best.pt / last.pt: checkpoint with the joint state dict, both configs, step and metrics;
* recon_val.png / recon_train.png: true and decoded CDF and DCDF per layout;
* latent_stats.json: SQ-Diffuse's target_mu / target_std (Sec. 3.4);
* log.txt.

Downstream code calls load_sq_vae(), encode_layout() and decode_fields().
"""

import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np
import torch
import torch.nn as nn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geom_ae_2d as gae                                              # noqa: E402
from geom_ae_2d import Attention, FeedForward, PointEmbed, PreNorm   # noqa: E402

REPO = gae.REPO
DEFAULT_DATA = os.path.join(REPO, "py", "squadgen", "data", "cdf")
DEFAULT_GEOM = os.path.join(REPO, "py", "squadgen", "runs", "geom_ae_2d", "best.pt")
# The manifest's method names: the five, and "oracle" (cdf_data.py: ORACLE's choice among them).
METHODS = ("zipline", "umber", "meridian", "torsion", "atlas", "oracle")
# The seven offsets in the order of the NPZ's corner ids: offset_abcd_id (own chart corner =
# dual-chart centre, side mid-point, chart centre, side mid-point), then offset_123_id (the
# chart's other three corners). Sec. 2.4.
OFFSET_KEYS = ("offset", "offsetb", "offsetc", "offsetd", "offset1", "offset2", "offset3")
IN_KEYS = ("dcdf", "gdcdf", "cdf", "gcdf", "sdf") + OFFSET_KEYS
OUT_KEYS = ("dcdf", "cdf", "offset", "offset1", "offset2", "offset3", "offsetb", "offsetc", "offsetd")
EVAL_QUERIES = 4096
NPZ_KEYS = ("subdiv_vertex", "subdiv_facet", "quad_split", "face2quad", "offset_abcd_id", "offset_123_id",
            "edge_vertex_color", "edge_color_store", "sample_bary_coords", "sample_point_tri_id")


# =============================================================================
# 1. FIELDS  (numpy; a layout is CDFGen's split quad mesh, triangulated)
# =============================================================================


def cross2(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def perp(a):
    return np.stack([a[..., 1], -a[..., 0]], -1)


class Layout:
    """One NPZ's split mesh as triangles, with the vertices in the Geom-AE's frame. Each
    triangle carries what its fields need: the coordinate pairs of its two quad edges (Sec. 2.2
    step 5) and the vertex ids of the seven offset targets. `centre`/`scale` map the frame back,
    x_obj = x / scale + centre, as geom_ae_2d.Domain does."""

    def __init__(self, path, fps_sizes=()):
        self.path = path
        with np.load(path) as z:
            a = {k: z[k] for k in NPZ_KEYS}
            self.fps = {n: (z[f"fps_point_tri_id_{n}"][..., 0].astype(np.int64),
                            z[f"fps_bary_coords_{n}"].astype(np.float64)) for n in fps_sizes}
        Q = a["subdiv_facet"].astype(np.int64)
        split = a["quad_split"].ravel().astype(bool)
        T = Q[:, [0, 1, 2, 2, 3, 0]]
        T[split] = Q[split][:, [3, 0, 1, 1, 2, 3]]
        self.T = T.reshape(-1, 3)
        # A triangle's edges v0->v1 and v1->v2 are quad edges, the third is the diagonal. Each
        # quad half-edge holds a signed index into the (start, end) pairs; a negative index is
        # the pair read backwards (load_patch_data_reformat, Appendix B).
        e = np.tile([0, 1, 2, 3], (len(Q), 1))
        e[split] = [3, 0, 1, 2]
        e = e.reshape(-1, 2) + 4 * (np.arange(len(self.T)) // 2)[:, None]
        cid = a["edge_vertex_color"].ravel()[e]
        col = a["edge_color_store"].astype(np.float64)[np.abs(cid)]
        col[cid < 0] = col[cid < 0][:, ::-1]
        self.col = col
        sub = a["face2quad"].ravel()[np.arange(len(self.T)) // 2]
        self.corners = np.concatenate([a["offset_abcd_id"], a["offset_123_id"]], 1).astype(np.int64)[sub]
        self.E = gae.boundary_edges(self.T)
        self.bverts = np.unique(self.E)
        V = a["subdiv_vertex"][:, :2].astype(np.float64)
        self.centre, self.scale = gae.normalisation(V, self.T, self.bverts)
        self.V = (V - self.centre) * self.scale
        self.samples = (a["sample_point_tri_id"].ravel().astype(np.int64), a["sample_bary_coords"].astype(np.float64))


def layout_fields(L, V, tri, bary):
    """Every point feature SQ-VAE reads, at the points (tri, bary) of layout L with vertices V.
    dcdf, gdcdf, cdf and gcdf are Appendix B's: the release's get_patch_grading_colors with
    scalar cross products, equal to it to float32 rounding. The two subchart coordinates c0, c1
    are oblique projections onto the triangle's two quad edges, mapped to their coordinate
    intervals; then CDF = 1 - max(c0, c1) and DCDF = min(c0, c1). The seven offsets are corner
    minus point. Reflections need no care: both cross products in t0, t1 change sign together,
    and the gradients turn with the mesh."""
    t = L.T[tri]
    v0, v1, v2 = V[t[:, 0]], V[t[:, 1]], V[t[:, 2]]
    p = bary[:, 0:1] * v0 + bary[:, 1:2] * v1 + bary[:, 2:3] * v2
    e0, e1 = v1 - v0, v2 - v1
    den = cross2(e0, e1)
    col = L.col[tri]
    t0 = np.clip(cross2(p - v0, e1) / den, 0, 1)            # along e0, parallel to e1
    t1 = np.clip(cross2(p - v1, e0) / -den, 0, 1)           # along e1, parallel to e0
    c0 = col[:, 0, 0] + t0 * (col[:, 0, 1] - col[:, 0, 0])
    c1 = col[:, 1, 0] + t1 * (col[:, 1, 1] - col[:, 1, 0])
    g0 = ((col[:, 0, 1] - col[:, 0, 0]) / den)[:, None] * perp(e1)
    g1 = ((col[:, 1, 1] - col[:, 1, 0]) / -den)[:, None] * perp(e0)
    g0 /= np.linalg.norm(g0, axis=1, keepdims=True) + 1e-10
    g1 /= np.linalg.norm(g1, axis=1, keepdims=True) + 1e-10
    out = {"xyz": p,
           "dcdf": np.clip(np.minimum(c0, c1), 0, 1)[:, None],
           "gdcdf": np.where((c1 < c0)[:, None], g1, g0),
           "cdf": np.clip(1 - np.maximum(c0, c1), 0, 1)[:, None],
           "gcdf": np.where((1 - c1 < 1 - c0)[:, None], -g1, -g0)}
    off = V[L.corners[tri]] - p[:, None]
    for i, k in enumerate(OFFSET_KEYS):
        out[k] = off[:, i]
    return out


# =============================================================================
# 2. DATA
# =============================================================================


def read_manifest(data, methods):
    path = os.path.join(data, "manifest.csv")
    if not os.path.exists(path):
        raise SystemExit(f"no {path}: run py/squadgen/make_cdf_data.sh first")
    with open(path) as f:
        rows = [r for r in csv.DictReader(f) if r["status"] == "ok" and r["method"] in methods]
    if not rows:
        raise SystemExit(f"no ok layouts of {' '.join(methods)} in {path}")
    for r in rows:
        r["npz"] = os.path.join(data, r["npz"])
    return rows


def make_sample(L, cfg, rng, aug):
    """One example under the release's batch keys (load_patch_data_reformat, 2D).
    Encoder input : xyz_fps_{N} and <k>_fps_{N} for every k in IN_KEYS (tokens), the same at
                    4N (keys/values), one random token set of each size out of the NPZ's five.
    Targets       : xyz_query and <k>_query for every k in OUT_KEYS, at n_query of the NPZ's
                    area samples."""
    V = gae.normalise(gae.augment(L.V, L.T, rng), L.T, L.bverts) if aug else L.V
    A, B = V[L.E[:, 0]], V[L.E[:, 1]]
    out = {}
    for n in (cfg["n_tokens"], 4 * cfg["n_tokens"]):
        tri, bary = L.fps[n]
        c = rng.integers(len(tri))
        f = layout_fields(L, V, tri[c], bary[c])
        f["sdf"] = gae.truncated_sdf(f["xyz"], A, B, cfg["tau"], known_inside=True)
        out.update({f"{k}_fps_{n}": v for k, v in f.items()})
    tri, bary = L.samples
    idx = rng.choice(len(tri), cfg["n_query"], replace=cfg["n_query"] > len(tri))
    f = layout_fields(L, V, tri[idx], bary[idx])
    out["xyz_query"] = f["xyz"]
    out.update({f"{k}_query": f[k] for k in OUT_KEYS})
    return {k: v.astype(np.float32) for k, v in out.items()}


class LayoutStream(torch.utils.data.IterableDataset):
    """Endless stream of augmented samples drawn uniformly over the training layouts. Each
    DataLoader worker seeds its own generator from (seed, worker id)."""

    def __init__(self, paths, cfg):
        self.paths, self.cfg = paths, cfg

    def __iter__(self):
        info = torch.utils.data.get_worker_info()
        rng = np.random.default_rng([self.cfg["seed"], info.id if info is not None else 0])
        sizes = (self.cfg["n_tokens"], 4 * self.cfg["n_tokens"])
        while True:
            L = Layout(self.paths[rng.integers(len(self.paths))], sizes)
            yield make_sample(L, self.cfg, rng, self.cfg["augment"])


# =============================================================================
# 3. NETWORK  (the release's KLAutoEncoder as SQ-VAE, models_ae.py, with the Sec. 3.1 edits)
# =============================================================================


class SQVAE(nn.Module):
    """The release's KLAutoEncoder with sqvae_geomae.yaml's latent_vae_params:
    is_context_norm = 0 (so no PreNorm normalises a context; the decoder still normalises its
    queries), use_geom_code_as_fps = use_geom_code_and_pe_as_fps = 1, gcolor_grad_norm_in = 1.
    No learnable latents, no encoder self-attention, no decoder FeedForward, no QK-norm.
    proj_geom_code_fps maps geom_dim to dim; the release has both 512, and here both are 256."""

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
            nn.Linear(dim, dim), nn.ReLU(), nn.Linear(dim, dim), nn.ReLU(), nn.Linear(dim, gae.KEY_DIM_DICT[k]))
            for k in self.out_keys_list_mlp})
        self.proj = nn.Linear(latent_dim, dim)
        self.mean_fc = nn.Linear(dim, latent_dim)
        self.logvar_fc = nn.Linear(dim, latent_dim)
        self.proj_geom_code_fps = nn.Linear(geom_dim, dim)
        self.initialize_weights()

    def initialize_weights(self):
        def basic(m):
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
        self.apply(basic)
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
        """-> {"latent": (B, N, latent_dim), "kl": (B,), "mean", "logvar"}. geom_code is the
        frozen Geom-AE's (B, N, geom_dim) code at xyz_fps_N (geom_code())."""
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
        """Every output key at queries (B, P, 2): what QuadExtraction will read off a
        triangulation of Omega (Sec. 4 step 6)."""
        return self.query_last_features(self.return_last_features(z), queries)


def geom_code(geom, batch):
    """The frozen Geom-AE's code at xyz_fps_N from its posterior mean (use_geom_cond_mean=1)."""
    with torch.no_grad():
        return geom.return_last_features_and_norm(geom.encode(batch, is_mode=True)["latent"])


# =============================================================================
# 4. TRAINING
# =============================================================================


MODEL_KEYS = ("depth", "dim", "geom_dim", "num_latents", "latent_dim", "heads", "dim_head", "drop_path_rate")


def build_model(cfg):
    return SQVAE(**{k: cfg[k] for k in MODEL_KEYS})


def calc_loss(model, geom, b, kl_weight, is_mode=False):
    """engine_ae.calc_loss with a frozen Geom-AE: mean L1 per output key (per sample over
    points and channels, then over the batch), weighted, plus kl_weight x KL."""
    enc = model.encode(b, geom_code(geom, b), is_mode=is_mode)
    pred = model.decode(enc["latent"], b["xyz_query"])
    terms = {k: (pred[k] - b[f"{k}_query"]).abs().mean(dim=[1, 2]).mean() for k in model.out_keys_list}
    terms["kl"] = enc["kl"].mean()
    loss = sum(w * terms[k] for k, w in zip(model.out_keys_list, model.out_weights_list))
    return loss + kl_weight * terms["kl"], terms


class LayoutEval:
    """Fixed, unaugmented samples of a layout set, decoded from the posterior mean at
    EVAL_QUERIES area samples each. Per layout:
      l1        weighted mean over the output keys of their L1 (best.pt is chosen on its mean)
      l1_<key>  each key's L1: dcdf and cdf are in [0, 1], offsets in the frame's units
      wall_iou  IoU of the decoded walls {cdf <= wall} with the true ones. wall = 0.1 is where
                QuadExtraction's region growing stops (-t 0.1, Sec. 5.1 step 2), so this is
                how well the decoded chart sides sit where the true ones are."""

    def __init__(self, rows, cfg):
        self.rows, self.items, self.wall = rows, [], cfg["wall"]
        small = dict(cfg, n_query=EVAL_QUERIES)
        for i, r in enumerate(rows):
            L = Layout(r["npz"], (cfg["n_tokens"], 4 * cfg["n_tokens"]))
            self.items.append(make_sample(L, small, np.random.default_rng(1000 + i), False))

    @torch.no_grad()
    def __call__(self, model, geom, dev, keep=False):
        was = model.training
        model.eval()
        per, latents = [], []
        w = dict(zip(model.out_keys_list, model.out_weights_list))
        for s in self.items:
            b = gae.to_dev({k: v[None] for k, v in s.items()}, dev)
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


def chart_sides(path, V):
    """Segments (S, 2, 2) of the true chart sides: split-mesh edges on dOmega or between
    faces of two charts (subchart_block, which BlockDecomposition.cdf() writes)."""
    with np.load(path) as z:
        Q = z["subdiv_facet"].astype(np.int64)
        chart = z["subchart_block"][z["face2quad"].ravel()]
    key = np.sort(np.concatenate([Q[:, [0, 1]], Q[:, [1, 2]], Q[:, [2, 3]], Q[:, [3, 0]]]), 1)
    face_chart = np.tile(chart, 4)
    uniq, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    lo = np.full(len(uniq), np.iinfo(np.int64).max)
    hi = np.full(len(uniq), np.iinfo(np.int64).min)
    np.minimum.at(lo, inv, face_chart)
    np.maximum.at(hi, inv, face_chart)
    return V[uniq[(cnt == 1) | (lo != hi)]]


def plot_reconstructions(ev, per, latents, model, dev, out, title, k=6):
    """A row per layout, worst L1 first, so that the failures are the first thing seen; one
    layout per domain, since the methods often agree on a domain's layout. The columns are the
    true and decoded CDF and DCDF, coloured per triangle of the split mesh at its centroid,
    with the true chart sides drawn over each. A row is as tall as its domain's aspect."""
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
def latent_stats(model, geom, paths, cfg, dev, n):
    """SQ-Diffuse's target_mu / target_std (Sec. 3.4): the per-channel mean and std of the
    posterior mean over n samples of the training layouts, augmented as SQ-Diffuse will see
    them. The release's values are for its 3D SQ-VAE and must not be reused."""
    model.eval()
    rng = np.random.default_rng(cfg["seed"] + 1)
    sizes = (cfg["n_tokens"], 4 * cfg["n_tokens"])
    s1 = s2 = 0.0
    count = 0
    for _ in range(math.ceil(n / cfg["batch"])):
        samples = [make_sample(Layout(paths[rng.integers(len(paths))], sizes), dict(cfg, n_query=1), rng,
                               cfg["augment"]) for _ in range(cfg["batch"])]
        b = gae.to_dev({k: np.stack([s[k] for s in samples]) for k in samples[0]}, dev)
        z = model.encode(b, geom_code(geom, b), is_mode=True)["latent"].cpu().double().flatten(0, 1)
        s1, s2, count = s1 + z.sum(0), s2 + (z * z).sum(0), count + len(z)
    mu = s1 / count
    std = (s2 / count - mu * mu).clamp_min(0).sqrt()
    return {"target_mu": mu.cpu().tolist(), "target_std": std.cpu().tolist(), "samples": n, "tokens": count}


def save(path, model, geom, cfg, gcfg, step, metrics):
    """main_ae.py's layout: KLAutoEncoderJoint's state dict, the Geom-AE under `cond_vae.` and
    SQ-VAE under `latent_vae.`."""
    sd = {f"cond_vae.{k}": v.detach().cpu() for k, v in geom.state_dict().items()}
    sd.update({f"latent_vae.{k}": v.detach().cpu() for k, v in model.state_dict().items()})
    torch.save({"model": sd, "config": gae.portable(cfg), "geom_config": gae.portable(gcfg), "step": step,
                "metrics": metrics}, path)


def load_sq_vae(path, device="cpu"):
    """SQ-VAE and its frozen Geom-AE from one checkpoint: returns (model, geom, config)."""
    ck = torch.load(path, map_location=device, weights_only=False)
    geom = gae.build_model(ck["geom_config"]).to(device)
    geom.load_state_dict({k[len("cond_vae."):]: v for k, v in ck["model"].items() if k.startswith("cond_vae.")})
    geom.requires_grad_(False)
    model = build_model(ck["config"]).to(device)
    model.load_state_dict({k[len("latent_vae."):]: v for k, v in ck["model"].items() if k.startswith("latent_vae.")})
    return model.eval(), geom.eval(), ck["config"]


@torch.no_grad()
def encode_layout(model, geom, cfg, npz_path, device="cpu", seed=0):
    """The posterior-mean latent (1, N, latent_dim) of one layout NPZ, unaugmented, in the
    Geom-AE's frame. Returns (latent, batch, Layout); Layout.centre/scale map the frame back."""
    L = Layout(npz_path, (cfg["n_tokens"], 4 * cfg["n_tokens"]))
    s = make_sample(L, dict(cfg, n_query=1), np.random.default_rng(seed), False)
    b = gae.to_dev({k: v[None] for k, v in s.items()}, device)
    return model.encode(b, geom_code(geom, b), is_mode=True)["latent"], b, L


@torch.no_grad()
def decode_fields(model, z, xy, chunk=8192):
    """Decode a latent (1, N, latent_dim) at points xy (P, 2) of the normalised frame, in
    chunks: {key: (P, channels) numpy}."""
    dev = z.device
    parts = [model.decode(z, torch.as_tensor(xy[s:s + chunk], dtype=torch.float32, device=dev)[None])
             for s in range(0, len(xy), chunk)]
    return {k: torch.cat([p[k][0] for p in parts]).cpu().numpy() for k in parts[0]}


def train(cfg, gcfg, geom, train_paths, val_eval, train_eval, log):
    dev = cfg["device"]
    amp = gae.amp_dtype(dev, cfg["amp"])
    torch.manual_seed(cfg["seed"])
    model = build_model(cfg).to(dev)
    geom = geom.to(dev)
    log(f"device={dev}  amp={amp}  workers={cfg['workers']}  SQ-VAE parameters="
        f"{sum(p.numel() for p in model.parameters()):,}  (Geom-AE frozen: {sum(p.numel() for p in geom.parameters()):,})")
    opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    scaler = gae.grad_scaler(dev, amp)
    loader = gae.make_loader(LayoutStream(train_paths, cfg), cfg)
    best, t0, run = float("inf"), time.time(), {}
    for step, batch in enumerate(loader, 1):
        for g in opt.param_groups:
            g["lr"] = gae.lr_at(step - 1, cfg)
        # The frozen Geom-AE runs under the same autocast, as the release runs both halves.
        with gae.autocast(dev, amp):
            loss, terms = calc_loss(model, geom, gae.to_dev(batch, dev), cfg["kl_weight"])
        gnorm = gae.train_step(model, opt, scaler, loss, cfg["clip_grad"])
        # The log line keeps the two fields apart and folds the seven offsets into one mean.
        shown = {"loss": loss, "dcdf": terms["dcdf"], "cdf": terms["cdf"],
                 "offsets": sum(terms[k] for k in OFFSET_KEYS) / len(OFFSET_KEYS), "kl": terms["kl"], "gnorm": gnorm}
        gae.accumulate(run, shown)
        if step % cfg["log_every"] == 0 or step == 1:
            el = time.time() - t0
            log(f"step {step:6d}  lr {opt.param_groups[0]['lr']:.2e}  {gae.averages(run)}"
                + f"  ({el / step * 1000:.0f} ms/step, {el / 60:.1f} min{gae.memory_note(dev)})")
            run = {}
        if step % cfg["eval_every"] == 0 or step == cfg["steps"]:
            mv, mt = val_eval(model, geom, dev), train_eval(model, geom, dev)
            log(f"  eval {step:6d}  val l1 {mv['l1']:.4f} cdf {mv['l1_cdf']:.4f} dcdf {mv['l1_dcdf']:.4f} "
                f"wall_iou {mv['wall_iou']:.4f}  |  train l1 {mt['l1']:.4f} cdf {mt['l1_cdf']:.4f} "
                f"dcdf {mt['l1_dcdf']:.4f} wall_iou {mt['wall_iou']:.4f}")
            save(os.path.join(cfg["out"], "last.pt"), model, geom, cfg, gcfg, step, {"val": mv, "train": mt})
            if mv["l1"] < best:
                best = mv["l1"]
                save(os.path.join(cfg["out"], "best.pt"), model, geom, cfg, gcfg, step, {"val": mv, "train": mt})
            gae.empty_cache(dev)
        if step >= cfg["steps"]:
            break


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--data", default=DEFAULT_DATA, help="make_cdf_data.sh's output (manifest.csv + npz/)")
    ap.add_argument("--geom_ae", default=DEFAULT_GEOM, help="geom_ae_2d.py checkpoint, frozen")
    ap.add_argument("--out", default=os.path.join(REPO, "py", "squadgen", "runs", "sq_vae_2d"))
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
    ap.add_argument("--workers", type=int, help="DataLoader workers (default: 2 on MPS/CPU, up to 8 on CUDA)")
    ap.add_argument("--device", default="auto", help="cuda, cuda:N, mps or cpu (default: the first there is)")
    ap.add_argument("--amp", choices=("auto", "bf16", "fp16", "none"), default="auto",
                    help="training autocast (auto: bf16 on CUDA that has it, else fp32)")
    ap.add_argument("--log_every", type=int, default=100)
    ap.add_argument("--eval_every", type=int, default=1000)
    ap.add_argument("--stats_samples", type=int, default=1024, help="samples for latent_stats.json")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--eval_only", action="store_true", help="evaluate, plot and write latent stats of best.pt")
    args = ap.parse_args()

    if not os.path.exists(args.geom_ae):
        raise SystemExit(f"no Geom-AE checkpoint {args.geom_ae}: train py/squadgen/geom_ae_2d.py first")
    geom, gcfg = gae.load_geom_ae(args.geom_ae)
    # Token i of SQ-VAE reads the point Geom-AE token i was computed at, so N is the Geom-AE's
    # (Sec. 3.2, "Token alignment"), and sdf is truncated at its tau.
    dev = gae.pick_device(args.device)
    cfg = dict(vars(args), n_tokens=gcfg["num_latents"], num_latents=gcfg["num_latents"], geom_dim=gcfg["dim"],
               tau=gcfg["tau"], augment=not args.no_aug, device=dev,
               workers=gae.default_workers(dev) if args.workers is None else args.workers)

    os.makedirs(args.out, exist_ok=True)
    logf = open(os.path.join(args.out, "log.txt"), "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  " + " ".join(f"{k}={v}" for k, v in sorted(cfg.items())))
    rows = read_manifest(args.data, args.methods)
    train_rows = [r for r in rows if r["split"] == "train"]
    val_rows = [r for r in rows if r["split"] == "val"]
    log(f"{len(rows)} layouts of {len({r['model'] for r in rows})} domains: {len(train_rows)} train, "
        f"{len(val_rows)} val ({' '.join(sorted({r['family'] for r in val_rows}))})")
    rng = np.random.default_rng(args.seed)
    val_eval = LayoutEval(val_rows, cfg)
    train_eval = LayoutEval([train_rows[i] for i in rng.choice(len(train_rows), min(48, len(train_rows)), replace=False)],
                            cfg)
    train_paths = [r["npz"] for r in train_rows]

    if not args.eval_only:
        train(cfg, gcfg, geom, train_paths, val_eval, train_eval, log)
    best = os.path.join(args.out, "best.pt")
    model, geom, _ = load_sq_vae(best, dev)
    step = torch.load(best, map_location="cpu", weights_only=False)["step"]
    for name, ev in (("val", val_eval), ("train", train_eval)):
        m, per, latents = ev(model, geom, dev, keep=True)
        log(f"best.pt (step {step}) on {name}: l1 {m['l1']:.4f}  wall_iou {m['wall_iou']:.4f}  "
            + "  ".join(f"{k} {m[f'l1_{k}']:.4f}" for k in OUT_KEYS)
            + f"  (wall IoU min {np.min([p['wall_iou'] for p in per]):.3f}, "
              f"median {np.median([p['wall_iou'] for p in per]):.3f})")
        plot_reconstructions(ev, per, latents, model, dev, os.path.join(args.out, f"recon_{name}.png"),
                             f"SQ-VAE 2D, best.pt (step {step}), {name} layouts: worst and best L1")
    stats = latent_stats(model, geom, train_paths, cfg, dev, args.stats_samples)
    with open(os.path.join(args.out, "latent_stats.json"), "w") as f:
        json.dump(dict(stats, step=step), f, indent=1)
    log(f"latent stats over {stats['tokens']} tokens: |mu| max {np.abs(stats['target_mu']).max():.3f}, "
        f"std {np.min(stats['target_std']):.3f}..{np.max(stats['target_std']):.3f} -> latent_stats.json")
    log(f"plots in {args.out}")


if __name__ == "__main__":
    main()
