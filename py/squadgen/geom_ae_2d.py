"""
geom_ae_2d.py -- train SQuadGen's Geometry Autoencoder (Geom-AE) on planar domains.

This is the 2D Geom-AE of docs/squadgen_2d.md Sec. 3.2, trained on the CrossGen corpus
(data/meshes/singlemat and data/meshes/mambo: 442 single-material domains) instead of the
random polygons of geom_ae_2d_prototype.py.

The network is a port of the release's KLAutoEncoder (SQuadGen squadgen/network/
models_ae.py) with the Sec. 3.1 edits for 2D: positions are 2-d (POS_DIM = 2) and the one
point feature is the truncated SDF where the release has the normal. Every module keeps
the release's name and shape, and the state dict is saved with the `cond_vae.` prefix that
main_ae.py writes, so a checkpoint loads into the edited release class as the frozen
geometry encoder of SQ-VAE (train_sqvae.sh --pretrained_geom_vae_pth, which keeps the keys
that start with `cond_vae.`). What is kept from the release, and why:

* Encoder: N token points (xyz_fps_N) cross-attend to 4N key/value points (xyz_fps_4N)
  through one head as wide as the model, then a FeedForward, then a diagonal Gaussian
  bottleneck of `latent_dim` channels. The bottleneck is sampled in training with KL
  weight 0 (Sec. 1.3, "Geom-AE"), so it is a noisy narrow layer rather than a prior.
* proj back to `dim`, `depth` pre-norm self-attention blocks (8 heads x 64 whatever `dim`
  is, because create_vae_from_config hard-codes it; Sec. 3.5), and the decoder's context
  LayerNorm. Its output, `return_last_features_and_norm`, is the N x dim geometry code that
  SQ-VAE and SQ-Diffuse consume frozen. Token i sits at xyz_fps_N[i], the point SQ-VAE's
  query i reads (Sec. 3.2, "Token alignment").
* Decoder: one single-head cross-attention from the query point's Fourier embedding to the
  tokens, then a 3-layer MLP head per output key, zero-initialised.

What changes for 2D (Sec. 3.2, "2D", and the Appendix D loader):

* Tokens and keys/values are AREA samples of Omega, not boundary samples. SQ-VAE's query
  points are area samples too (the NPZ's fps sets lie on the split quad mesh, Sec. 2.4),
  and token i has to sit at query point i.
* The feature is sdf = clip(phi / tau, -1, 1), tau = 0.1, phi < 0 inside (Sec. 2.7 item 3).
  On area samples it lies in [-1, 0] and tells a token how far, up to tau, it is from
  dOmega. Domain shape is encoded both by where the points are and by these values.
* The target is the same truncated SDF at three groups of query points: output keys
  `sdf` (area samples, the release's "surface" group `xyz_query`), `sdf_near` (boundary
  samples + N(0, 0.02^2)) and `sdf_global` (uniform in [-1, 1]^2), weighted (1, 1, 0.1)
  (Sec. 3.1 item 7). The truncation is what keeps the 0.1 box weight from collapsing the
  output to phi = 0, which an untruncated SDF did (Sec. 3.2, "Collapse risk"); the release
  truncates its UDF at 0.08 for the same reason. --w_global 1 is the fallback.

Data. Like the release's NPZs (Sec. 2.4, `fps_bary_coords_{n}`), each domain stores five
farthest-point-sampled token sets per size as (triangle, barycentric) pairs. They are
computed once, in the normalised frame, and cached. Because they are barycentric they
follow the mesh through any deformation. Augmentation is the 2D one of Sec. 2.7 item 1:
rotation about z, reflection, anisotropic scale and a smooth warp. A deformation that folds
a triangle is rejected. The inference normalisation of Sec. 4 step 1 follows (centre, then
0.999 / max radius), so training sees domains in the frame inference will hand it. Mambo
names are <model>_<face> (S39_12 is planar face 12 of CAD model S39). The train/val split
is by model, so validation measures generalisation to shapes from a model the network has
not seen in any face.

The release trainers need CUDA and torchrun (Sec. 3.4, "Hardware"), so this is the plain
loop around calc_loss that the guide suggests. It runs on CUDA, Apple MPS or the CPU.

Run
---
    python py/squadgen/geom_ae_2d.py                                  # train -> py/squadgen/runs/geom_ae_2d
    python py/squadgen/geom_ae_2d.py --out /tmp/smoke --steps 200     # smoke test
    python py/squadgen/geom_ae_2d.py --eval_only                      # plots from best.pt

Memory: on an 8 GB Apple M2, batch 8 with the default query counts takes ~3.3 GB of MPS
memory at ~0.5 s/step; batch 12 needs 4.3 GB and batch 16 swapped (100 s/step). bf16
autocast was slower on MPS, so there is no AMP option.

Outputs in --out: fps_cache_N{N}.npz (token sets), best.pt / last.pt (checkpoint: model
under `cond_vae.`, config, step, metrics), recon_val.png / recon_train.png (decoded SDF
against the true boundary), log.txt. Downstream code calls load_geom_ae() and
encode_obj().
"""

import argparse
import glob
import math
import os
import time
import zlib

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_DIRS = [os.path.join(REPO, "data", "meshes", "singlemat"),
                os.path.join(REPO, "data", "meshes", "mambo")]
N_FPS_COPIES = 5        # token sets per domain and size, as the release's NPZs store
POS_DIM = 2             # Sec. 3.1 item 1
# Channels of every point feature and output key: the release's network/utils.py KEY_DIM_DICT
# with the Sec. 3.1 item 1 edits. One table for both networks, as there: the Geom-AE reads the
# sdf keys, SQ-VAE (sq_vae_2d.py) the fields, their unit gradients and the seven offsets,
# every vector of them 2-d.
KEY_DIM_DICT = {"sdf": 1, "sdf_near": 1, "sdf_global": 1, "dcdf": 1, "cdf": 1, "gdcdf": 2, "gcdf": 2,
                "offset": 2, "offset1": 2, "offset2": 2, "offset3": 2, "offsetb": 2, "offsetc": 2, "offsetd": 2}
GEOM_KEYS = ["sdf"]


# =============================================================================
# 1. GEOMETRY  (numpy; a domain is a triangle mesh of Omega at z = 0)
# =============================================================================


def read_obj(path):
    """Vertices (V, 2) and triangles (F, 3) of a CrossGen OBJ. Material tags are ignored:
    the Geom-AE sees only the domain."""
    V, Fc = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                V.append(line.split()[1:3])
            elif line.startswith("f "):
                Fc.append([int(t.split("/")[0]) - 1 for t in line.split()[1:4]])
    return np.asarray(V, np.float64), np.asarray(Fc, np.int64)


def signed_areas(V, Fc):
    a, b, c = V[Fc[:, 0]], V[Fc[:, 1]], V[Fc[:, 2]]
    return 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))


def boundary_edges(Fc):
    """Edges used by exactly one triangle. They close up into the loops of dOmega, so the
    even-odd test below needs no orientation."""
    E = np.concatenate([Fc[:, [0, 1]], Fc[:, [1, 2]], Fc[:, [2, 0]]])
    _, inv, cnt = np.unique(np.sort(E, 1), axis=0, return_inverse=True, return_counts=True)
    return E[cnt[inv.ravel()] == 1]


def normalisation(V, Fc, bverts):
    """The inference frame (Sec. 4 step 1) as (centre, scale): area centroid at the origin,
    farthest boundary vertex at radius 0.999. The area centroid rather than the vertex mean,
    so that the frame does not depend on how the domain was meshed."""
    ar = np.abs(signed_areas(V, Fc))
    c = (ar[:, None] * V[Fc].mean(1)).sum(0) / ar.sum()
    return c, 0.999 / np.linalg.norm(V[bverts] - c, axis=1).max()


def normalise(V, Fc, bverts):
    c, s = normalisation(V, Fc, bverts)
    return (V - c) * s


def point_from_bary(V, Fc, tri, bary):
    return np.einsum("pk,pkd->pd", bary, V[Fc[tri]])


def sample_area(V, Fc, cdf, n, rng):
    """n uniform samples of Omega as (triangle, barycentric) pairs. `cdf` is the cumulative
    triangle area; it is a ratio, so it survives the similarity part of any augmentation
    and the anisotropic part only changes the density by up to 10 %."""
    tri = np.minimum(np.searchsorted(cdf, rng.random(n) * cdf[-1]), len(cdf) - 1)
    u, v = rng.random(n), rng.random(n)
    flip = u + v > 1
    u[flip], v[flip] = 1 - u[flip], 1 - v[flip]
    return tri, np.stack([1 - u - v, u, v], 1)


def farthest_point_sampling(pts, n, rng):
    """Greedy max-min selection (the release uses torch_cluster.fps)."""
    chosen = np.empty(n, np.int64)
    chosen[0] = rng.integers(len(pts))
    d = np.full(len(pts), np.inf)
    for i in range(1, n):
        d = np.minimum(d, ((pts - pts[chosen[i - 1]]) ** 2).sum(1))
        chosen[i] = d.argmax()
    return chosen


def _chunk(n_edges):
    """Points per chunk so that a (chunk, E, 2) temporary stays near 16 MB. Several
    DataLoader workers each hold a few of them, on machines with 8 GB."""
    return max(16, (1 << 20) // max(1, n_edges))


def _point_segment(p, a, b):
    e = b - a
    t = np.clip(((p - a) * e).sum(-1) / np.maximum((e * e).sum(-1), 1e-30), 0.0, 1.0)
    return np.sqrt(((p - a - t[..., None] * e) ** 2).sum(-1))


def segment_distance(x, A, B, cap=np.inf):
    """Distance from each x (P, 2) to the segment set {[A_e, B_e]}: exact where it is below
    `cap`, and `cap` elsewhere. The truncated SDF needs nothing beyond tau, and a finite cap
    lets each point look only at the edges registered in its cell of a cap-sized grid, about
    a tenth of dOmega. That is what keeps two DataLoader workers ahead of the GPU. An edge
    is registered in every cell that its bounding box, grown by cap, touches. A point within
    cap of the edge lies in that grown box, so its cell has the edge."""
    if not np.isfinite(cap):
        out = np.empty(len(x))
        step = _chunk(len(A))
        for s in range(0, len(x), step):
            out[s:s + step] = _point_segment(x[s:s + step, None], A[None], B[None]).min(1)
        return out
    lo = np.minimum(x.min(0), np.minimum(A, B).min(0)) - cap
    emin = np.floor((np.minimum(A, B) - cap - lo) / cap).astype(np.int64)
    emax = np.floor((np.maximum(A, B) + cap - lo) / cap).astype(np.int64)
    pc = np.floor((x - lo) / cap).astype(np.int64)
    ny = max(emax[:, 1].max(), pc[:, 1].max()) + 1
    span = emax - emin + 1
    per = span[:, 0] * span[:, 1]
    eid = np.repeat(np.arange(len(A)), per)
    local = np.arange(per.sum()) - np.repeat(np.cumsum(per) - per, per)
    cell = (emin[eid, 0] + local % span[eid, 0]) * ny + emin[eid, 1] + local // span[eid, 0]
    order = np.argsort(cell, kind="stable")
    cell, eid = cell[order], eid[order]
    pcell = pc[:, 0] * ny + pc[:, 1]
    start = np.searchsorted(cell, pcell, "left")
    n = np.searchsorted(cell, pcell, "right") - start
    out = np.full(len(x), cap)
    if n.sum() == 0:
        return out
    pid = np.repeat(np.arange(len(x)), n)
    e = eid[start[pid] + np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)]
    d = _point_segment(x[pid], A[e], B[e])
    has = n > 0
    out[has] = np.minimum(cap, np.minimum.reduceat(d, (np.cumsum(n) - n)[has]))
    return out


def inside(x, A, B):
    """Even-odd rule, ray in +x. Holes need no special treatment."""
    out = np.empty(len(x), bool)
    ya, yb = A[None, :, 1], B[None, :, 1]
    dy = np.where(yb == ya, 1.0, yb - ya)       # such edges never straddle; avoids 0/0
    step = _chunk(len(A))
    for s in range(0, len(x), step):
        y = x[s:s + step, 1:2]
        straddles = (ya > y) != (yb > y)
        xc = A[None, :, 0] + (y - ya) / dy * (B - A)[None, :, 0]
        out[s:s + step] = (straddles & (x[s:s + step, 0:1] < xc)).sum(1) % 2 == 1
    return out


def truncated_sdf(x, A, B, tau, known_inside=False):
    """clip(phi / tau, -1, 1) as a (P, 1) column, phi < 0 inside Omega. Area samples are
    inside by construction, and skipping the ray cast for them keeps a sample that lands on
    dOmega from being misclassified by rounding."""
    d = segment_distance(x, A, B, cap=tau)
    d = -d if known_inside else np.where(inside(x, A, B), -d, d)
    return np.clip(d / tau, -1.0, 1.0).astype(np.float32)[:, None]


def smooth_warp(p, rng, amp):
    """A few random low-frequency sinusoids. They play the part of the release's Perlin
    noise without needing the `noise` package."""
    out = p.copy()
    for _ in range(3):
        k = rng.normal(0, 1.5, 2)
        ph = rng.uniform(0, 2 * np.pi)
        out += amp * rng.normal(0, 1, 2) * np.sin(p @ k + ph)[:, None] / 3
    return out


def augment(V, Fc, rng, warp_amp=0.06):
    """Sec. 2.7 item 1 / Appendix D add_noise_2d, without the translation: normalise()
    re-centres anyway. A deformation that folds a triangle is rejected and retried. A
    reflection flips every triangle, which the test accepts."""
    a0 = signed_areas(V, Fc)
    for _ in range(5):
        p = V.copy()
        if rng.random() < 0.5:
            p = smooth_warp(p, rng, rng.uniform(0, warp_amp))
        th = rng.uniform(-np.pi, np.pi)
        p = p @ np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]).T
        if rng.random() < 0.5:
            p[:, 0] = -p[:, 0]
        p = p * rng.uniform(0.9, 1.1, 2)
        s = signed_areas(p, Fc) * a0
        if np.all(s > 0) or np.all(s < 0):
            return p
    return V.copy()


# =============================================================================
# 2. DATA
# =============================================================================


def family(path):
    """Mambo S39_12 -> 'S39'; singlemat geom012 -> 'geom012' (one family per model)."""
    name = os.path.splitext(os.path.basename(path))[0]
    return name.rsplit("_", 1)[0] if "_" in name else name


class Domain:
    """One OBJ in the normalised frame. `centre`/`scale` map it back: x_obj = x / scale +
    centre, the transform Sec. 4 step 1 says to keep."""

    def __init__(self, path=None, V=None, Fc=None):
        self.path = path
        self.name = os.path.basename(path) if path else "domain"
        if path is not None:
            V, Fc = read_obj(path)
        self.F = Fc.astype(np.int32)       # every DataLoader worker holds a copy of the corpus
        self.E = boundary_edges(self.F)
        self.bverts = np.unique(self.E)
        self.centre, self.scale = normalisation(V, Fc, self.bverts)
        self.V = (V - self.centre) * self.scale
        self.cdf = np.cumsum(np.abs(signed_areas(self.V, self.F)))
        self.fps = None      # {n: ((copies, n) triangle ids, (copies, n, 3) barycentrics)}


def fps_sets(d, sizes, rng):
    """N_FPS_COPIES token sets of each size for one domain: FPS on 8n (at least 8192) area
    samples. The release draws them from 50 000 Poisson-disk samples (Sec. 2.4)."""
    out = {}
    for n in sizes:
        tris, barys = [], []
        for _ in range(N_FPS_COPIES):
            tri, bary = sample_area(d.V, d.F, d.cdf, max(8 * n, 8192), rng)
            sel = farthest_point_sampling(point_from_bary(d.V, d.F, tri, bary), n, rng)
            tris.append(tri[sel])
            barys.append(bary[sel])
        out[n] = (np.stack(tris).astype(np.int32), np.stack(barys).astype(np.float32))
    return out


def _fps_job(job):
    path, sizes = job
    return path, fps_sets(Domain(path), sizes, np.random.default_rng(zlib.crc32(os.path.basename(path).encode())))


def load_domains(dirs, n_tokens, cache_dir, workers):
    """Load every OBJ and attach its cached FPS token sets, computing them on first use.
    A cache entry is keyed by file name and face count. A regenerated mesh with a different
    face count gets new sets, because its triangle ids would point elsewhere."""
    paths = sorted(p for d in dirs for p in glob.glob(os.path.join(d, "*.obj")))
    if not paths:
        raise SystemExit(f"no .obj files in {dirs}")
    names = [os.path.basename(p) for p in paths]
    if len(set(names)) != len(names):
        raise SystemExit("two domains share a file name; the FPS cache is keyed by it")
    domains = [Domain(p) for p in paths]
    sizes = (n_tokens, 4 * n_tokens)
    cache = os.path.join(cache_dir, f"fps_cache_N{n_tokens}.npz")
    stored = dict(np.load(cache)) if os.path.exists(cache) else {}
    todo = [d.path for d in domains
            if f"{d.name}/nf" not in stored or int(stored[f"{d.name}/nf"]) != len(d.F)]
    if todo:
        print(f"computing FPS token sets for {len(todo)} domains ({workers} workers) ...", flush=True)
        t0 = time.time()
        jobs = [(p, sizes) for p in todo]
        if workers > 1:
            import multiprocessing as mp
            with mp.get_context("spawn").Pool(workers) as pool:
                results = pool.map(_fps_job, jobs, chunksize=4)
        else:
            results = [_fps_job(j) for j in jobs]
        nf = {d.path: len(d.F) for d in domains}
        for p, out in results:
            name = os.path.basename(p)
            stored[f"{name}/nf"] = np.int64(nf[p])
            for n, (tri, bary) in out.items():
                stored[f"{name}/tri{n}"], stored[f"{name}/bary{n}"] = tri, bary
        os.makedirs(cache_dir, exist_ok=True)
        np.savez(cache, **stored)
        print(f"  done in {time.time() - t0:.0f} s -> {cache}", flush=True)
    for d in domains:
        d.fps = {n: (stored[f"{d.name}/tri{n}"], stored[f"{d.name}/bary{n}"]) for n in sizes}
    return domains


def split_by_family(domains, val_frac, seed):
    fams = sorted({family(d.path) for d in domains})
    rng = np.random.default_rng(seed)
    val = set(rng.choice(fams, max(1, round(val_frac * len(fams))), replace=False))
    return ([d for d in domains if family(d.path) not in val],
            [d for d in domains if family(d.path) in val])


def make_sample(d, cfg, rng, aug):
    """One example under the release's batch keys (Appendix D, Geom-AE part only).
    Encoder input : xyz_fps_{N}, sdf_fps_{N} (tokens); xyz_fps_{4N}, sdf_fps_{4N} (keys/values).
    Targets       : sdf_query at xyz_query (area samples), sdf_near at xyz_near, and
                    sdf_global at xyz_global, all truncated SDF in units of tau."""
    V = normalise(augment(d.V, d.F, rng), d.F, d.bverts) if aug else d.V
    A, B = V[d.E[:, 0]], V[d.E[:, 1]]
    tau, out = cfg["tau"], {}
    for n in (cfg["n_tokens"], 4 * cfg["n_tokens"]):
        tri, bary = d.fps[n]
        c = rng.integers(len(tri))
        x = point_from_bary(V, d.F, tri[c], bary[c].astype(np.float64))
        out[f"xyz_fps_{n}"] = x.astype(np.float32)
        out[f"sdf_fps_{n}"] = truncated_sdf(x, A, B, tau, known_inside=True)
    x = point_from_bary(V, d.F, *sample_area(V, d.F, d.cdf, cfg["n_query"], rng))
    out["xyz_query"], out["sdf_query"] = x.astype(np.float32), truncated_sdf(x, A, B, tau, True)
    L = np.linalg.norm(B - A, axis=1)
    ie = np.minimum(np.searchsorted(np.cumsum(L), rng.random(cfg["n_near"]) * L.sum()), len(L) - 1)
    x = (A[ie] + rng.random((len(ie), 1)) * (B[ie] - A[ie])
         + rng.normal(0, cfg["sigma_near"], (len(ie), 2)))
    out["xyz_near"], out["sdf_near"] = x.astype(np.float32), truncated_sdf(x, A, B, tau)
    x = rng.uniform(-1, 1, (cfg["n_global"], 2))
    out["xyz_global"], out["sdf_global"] = x.astype(np.float32), truncated_sdf(x, A, B, tau)
    return out


class DomainStream(torch.utils.data.IterableDataset):
    """Endless stream of augmented samples drawn uniformly over the training domains. Each
    DataLoader worker seeds its own generator from (seed, worker id)."""

    def __init__(self, domains, cfg):
        self.domains, self.cfg = domains, cfg

    def __iter__(self):
        info = torch.utils.data.get_worker_info()
        rng = np.random.default_rng([self.cfg["seed"], info.id if info is not None else 0])
        while True:
            d = self.domains[rng.integers(len(self.domains))]
            yield make_sample(d, self.cfg, rng, self.cfg["augment"])


# =============================================================================
# 3. NETWORK  (the release's KLAutoEncoder, models_ae.py, with the Sec. 3.1 edits)
# =============================================================================
# Module names and shapes are the release's, so that state_dict keys match it. The release
# also loops the key/value points over "groups" and max-pools the encoder output across
# them. A [B, 4N, 2] input is one group, which makes that a no-op, so it is left out.


class DropPath(nn.Module):
    """Per-sample stochastic depth (timm's DropPath, which the release imports)."""

    def __init__(self, p=0.0):
        super().__init__()
        self.p = p

    def forward(self, x):
        if self.p == 0.0 or not self.training:
            return x
        keep = 1.0 - self.p
        return x * x.new_empty((x.shape[0],) + (1,) * (x.ndim - 1)).bernoulli_(keep) / keep


class PreNorm(nn.Module):
    def __init__(self, dim, fn, context_dim=None, is_context_norm=1, is_x_norm=1):
        super().__init__()
        self.fn = fn
        self.norm = nn.LayerNorm(dim) if is_x_norm else nn.Identity()
        self.norm_context = nn.LayerNorm(context_dim) if context_dim is not None and is_context_norm else None
        self.is_context_norm = is_context_norm

    def forward(self, x, **kwargs):
        x = self.norm(x)
        if self.norm_context is not None:
            kwargs["context"] = self.norm_context(kwargs["context"])
        return self.fn(x, **kwargs)


class GEGLU(nn.Module):
    def forward(self, x):
        x, gates = x.chunk(2, dim=-1)
        return x * F.gelu(gates)


class FeedForward(nn.Module):
    def __init__(self, dim, mult=4, drop_path_rate=0.0):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(dim, dim * mult * 2), GEGLU(), nn.Linear(dim * mult, dim))
        self.drop_path = DropPath(drop_path_rate)

    def forward(self, x, **kwargs):
        return self.drop_path(self.net(x))


class Attention(nn.Module):
    def __init__(self, query_dim, context_dim=None, heads=8, dim_head=64, drop_path_rate=0.0):
        super().__init__()
        inner = heads * dim_head
        self.heads, self.scale = heads, dim_head ** -0.5
        self.to_q = nn.Linear(query_dim, inner, bias=False)
        self.to_kv = nn.Linear(context_dim or query_dim, 2 * inner, bias=False)
        self.to_out = nn.Linear(inner, query_dim)
        self.drop_path = DropPath(drop_path_rate)

    def forward(self, x, context=None):
        context = x if context is None else context
        k, v = self.to_kv(context).chunk(2, dim=-1)
        split = lambda t: t.unflatten(-1, (self.heads, -1)).transpose(1, 2)   # b n (h d) -> b h n d
        o = F.scaled_dot_product_attention(split(self.to_q(x)), split(k), split(v), scale=self.scale)
        return self.drop_path(self.to_out(o.transpose(1, 2).flatten(-2)))


class PointEmbed(nn.Module):
    """Learned Fourier features of the position. mlp_pos holds random frequencies
    (initialised N(0, 1), times 2 pi) and is trained. One linear layer per point feature,
    and a linear layer over the concatenation. `query` is the position part alone; the
    decoder embeds its query points with it, so they share mlp_pos and mlp_pos_post with
    the encoder."""

    def __init__(self, dim, in_keys_list=()):
        super().__init__()
        self.mlp_pos = nn.Linear(POS_DIM, dim // 2, bias=False)
        self.mlp_pos_post = nn.Linear(dim // 2 * 2 + POS_DIM, dim)
        for k in in_keys_list:
            setattr(self, f"mlp_{k}", nn.Sequential(nn.Linear(KEY_DIM_DICT[k], dim)))
        self.mlp = nn.Linear(dim * (1 + len(in_keys_list)), dim)

    def query(self, xyz):
        f = self.mlp_pos(xyz) * 2 * np.pi
        return self.mlp_pos_post(torch.cat([xyz, f.sin(), f.cos()], dim=-1))

    def forward(self, input, context_dict):
        h = [self.query(input)] + [getattr(self, f"mlp_{k}")(v) for k, v in context_dict.items()]
        return self.mlp(torch.cat(h, dim=-1))


class KLAutoEncoder(nn.Module):
    """The release's KLAutoEncoder as the Geom-AE uses it: is_context_norm = 1 (so the
    encoder cross-attention normalises its context and the decoder normalises both sides),
    no learnable latents, no encoder self-attention, no decoder FeedForward, no QK-norm."""

    def __init__(self, in_keys_list=("sdf",), out_keys_list=("sdf", "sdf_near", "sdf_global"),
                 out_weights_list=(1.0, 1.0, 0.1), depth=8, dim=256, num_latents=512,
                 latent_dim=32, heads=8, dim_head=64, drop_path_rate=0.1):
        super().__init__()
        self.num_latents = num_latents
        self.in_keys_list_mlp = sorted({k for k in in_keys_list})
        self.out_keys_list = list(out_keys_list)
        self.out_weights_list = list(out_weights_list)
        self.out_keys_list_mlp = sorted({k.replace("_near", "").replace("_global", "") for k in out_keys_list})

        self.cross_attend_blocks = nn.ModuleList([
            PreNorm(dim, Attention(dim, dim, heads=1, dim_head=dim), context_dim=dim),
            PreNorm(dim, FeedForward(dim)),
        ])
        self.point_embed = PointEmbed(dim, self.in_keys_list_mlp)
        self.layers = nn.ModuleList([nn.ModuleList([
            PreNorm(dim, Attention(dim, heads=heads, dim_head=dim_head, drop_path_rate=drop_path_rate)),
            PreNorm(dim, FeedForward(dim, drop_path_rate=drop_path_rate)),
        ]) for _ in range(depth)])
        self.decoder_cross_attn = PreNorm(dim, Attention(dim, dim, heads=1, dim_head=dim), context_dim=dim)
        self.to_outputs = nn.ModuleDict({k: nn.Sequential(
            nn.Linear(dim, dim), nn.ReLU(), nn.Linear(dim, dim), nn.ReLU(), nn.Linear(dim, KEY_DIM_DICT[k]))
            for k in self.out_keys_list_mlp})
        self.proj = nn.Linear(latent_dim, dim)
        self.mean_fc = nn.Linear(dim, latent_dim)
        self.logvar_fc = nn.Linear(dim, latent_dim)
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

    def encode(self, batch, is_mode=False):
        """-> {"latent": (B, N, latent_dim), "kl": (B,), "mean", "logvar"}. A sample of the
        posterior in training, its mean with is_mode (what SQ-VAE uses, use_geom_cond_mean=1)."""
        n = self.num_latents
        ctx = lambda s: {k: batch[f"{k}_fps_{s}"] for k in self.in_keys_list_mlp}
        x = self.point_embed(batch[f"xyz_fps_{n}"], ctx(n))
        pc = self.point_embed(batch[f"xyz_fps_{4 * n}"], ctx(4 * n))
        cross_attn, cross_ff = self.cross_attend_blocks
        x = cross_attn(x, context=pc) + x
        x = cross_ff(x) + x
        # DiagonalGaussianDistribution: logvar clamped to [-10, 10], KL averaged over
        # tokens and channels (dims [1, 2]).
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

    def return_last_features_and_norm(self, z):
        """The geometry code: what SQ-VAE's proj_geom_code_fps and SQ-Diffuse read."""
        return self.decoder_cross_attn.norm_context(self.return_last_features(z))

    def query_last_features(self, x, queries):
        """x from return_last_features (NOT the normed code: the decoder's PreNorm applies
        norm_context itself)."""
        h = self.decoder_cross_attn(self.point_embed.query(queries), context=x)
        return {k: self.to_outputs[k](h) for k in self.out_keys_list_mlp}

    def decode(self, z, queries):
        return self.query_last_features(self.return_last_features(z), queries)

    def decode_batch(self, z, batch):
        """Every output key at its own query group: sdf at xyz_query, sdf_near at xyz_near,
        sdf_global at xyz_global, decoded in one pass."""
        groups = [g for g in ("query", "near", "global") if any(group_of(k) == g for k in self.out_keys_list)]
        q = [batch[f"xyz_{g}"] for g in groups]
        out = self.decode(z, torch.cat(q, dim=1))
        bounds = np.cumsum([0] + [t.shape[1] for t in q])
        res = {}
        for k in self.out_keys_list:
            i = groups.index(group_of(k))
            res[k] = out[k.replace("_near", "").replace("_global", "")][:, bounds[i]:bounds[i + 1]]
        return res


def group_of(key):
    return "near" if key.endswith("_near") else "global" if key.endswith("_global") else "query"


def label_of(key):
    """The release's get_label: `<k>_near` / `<k>_global` are batch keys themselves, a
    surface key reads `<k>_query`."""
    return key if group_of(key) != "query" else f"{key}_query"


# =============================================================================
# 4. TRAINING
# =============================================================================


MODEL_KEYS = ("depth", "dim", "num_latents", "latent_dim", "heads", "dim_head", "drop_path_rate")


def build_model(cfg):
    return KLAutoEncoder(in_keys_list=GEOM_KEYS, out_keys_list=("sdf", "sdf_near", "sdf_global"),
                         out_weights_list=(cfg["w_query"], cfg["w_near"], cfg["w_global"]),
                         **{k: cfg[k] for k in MODEL_KEYS})


def pick_device():
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def to_dev(batch, dev):
    return {k: torch.as_tensor(v).to(dev, non_blocking=True) for k, v in batch.items()}


def calc_loss(model, b, kl_weight, is_mode=False):
    """engine_ae.calc_loss for the Geom-AE alone: mean L1 per output key (per sample over
    points and channels, then over the batch), weighted, plus kl_weight x KL."""
    enc = model.encode(b, is_mode=is_mode)
    pred = model.decode_batch(enc["latent"], b)
    terms = {k: (pred[k] - b[label_of(k)]).abs().mean(dim=[1, 2]).mean() for k in model.out_keys_list}
    terms["kl"] = enc["kl"].mean()
    loss = sum(w * terms[k] for k, w in zip(model.out_keys_list, model.out_weights_list))
    return loss + kl_weight * terms["kl"], terms


class GridEval:
    """Fixed, unaugmented evaluation of a domain set, using the posterior mean. The decoded
    SDF on a res^2 grid of [-1, 1]^2 is compared against the exact one, and three numbers
    are reported:
      band_l1  L1 error where the true truncated SDF is inside the band (|phi| < tau), in
               units of tau. That is where the target carries information.
      iou      IoU of the decoded interior {phi_hat < 0} with Omega.
      bnd      mean |phi_hat| on dOmega in units of tau: roughly how far the decoded zero
               level set sits from the true boundary, as a fraction of tau."""

    def __init__(self, domains, cfg, res=128):
        g = np.linspace(-1, 1, res)
        X, Y = np.meshgrid(g, g)
        self.grid = np.stack([X.ravel(), Y.ravel()], 1)
        self.res, self.items = res, []
        small = dict(cfg, n_query=16, n_near=16, n_global=16)
        for i, d in enumerate(domains):
            rng = np.random.default_rng(1000 + i)
            s = make_sample(d, small, rng, False)
            A, B = d.V[d.E[:, 0]], d.V[d.E[:, 1]]
            L = np.linalg.norm(B - A, axis=1)
            ie = np.minimum(np.searchsorted(np.cumsum(L), rng.random(2048) * L.sum()), len(L) - 1)
            on = A[ie] + rng.random((2048, 1)) * (B[ie] - A[ie])
            self.items.append((d, s, truncated_sdf(self.grid, A, B, cfg["tau"])[:, 0], on.astype(np.float32)))
        self.grid = self.grid.astype(np.float32)

    @torch.no_grad()
    def __call__(self, model, dev, keep=False):
        was = model.training
        model.eval()
        band, iou, bnd, preds = [], [], [], []
        for d, s, tgt, on in self.items:
            b = to_dev({k: v[None] for k, v in s.items()}, dev)
            x = model.return_last_features(model.encode(b, is_mode=True)["latent"])
            p = model.query_last_features(x, torch.from_numpy(self.grid)[None].to(dev))["sdf"][0, :, 0].cpu().numpy()
            pb = model.query_last_features(x, torch.from_numpy(on)[None].to(dev))["sdf"][0, :, 0].cpu().numpy()
            m = np.abs(tgt) < 1
            band.append(np.abs(p - tgt)[m].mean())
            iou.append(((p < 0) & (tgt < 0)).sum() / max(1, ((p < 0) | (tgt < 0)).sum()))
            bnd.append(np.abs(pb).mean())
            if keep:
                preds.append(p)
        model.train(was)
        res = {"band_l1": float(np.mean(band)), "iou": float(np.mean(iou)), "bnd": float(np.mean(bnd))}
        return (res, preds, iou) if keep else res


# Diverging field from sdf/tau = -1 to +1: red inside, neutral grey on dOmega, blue outside.
FIELD_CMAP = ("#e34948", "#f0efec", "#2a78d6")


def plot_reconstructions(ev, preds, ious, out, title, k=12):
    """A panel per domain, worst IoU first, so that the failures are the first thing seen.
    Each panel shows the decoded truncated SDF as filled colour, its zero level set as a
    solid line, the true dOmega dashed, and the tokens as grey dots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D

    order = np.argsort(ious)
    idx = np.concatenate([order[:k // 2], order[::-1][:k - k // 2]]) if len(order) > k else order
    cols = min(4, len(idx))
    rows = math.ceil(len(idx) / cols)
    cmap = LinearSegmentedColormap.from_list("sdf", FIELD_CMAP)
    fig, axes = plt.subplots(rows, cols, figsize=(3.3 * cols, 3.3 * rows + 0.6), squeeze=False)
    X = ev.grid[:, 0].reshape(ev.res, ev.res)
    Y = ev.grid[:, 1].reshape(ev.res, ev.res)
    cf = None
    for ax, i in zip(axes.ravel(), idx):
        d, s, _, _ = ev.items[i]
        P = preds[i].reshape(ev.res, ev.res)
        # The decoder is not clipped and overshoots +-1 a little where the target saturates;
        # contourf would leave those regions unfilled, which reads as an artefact.
        cf = ax.contourf(X, Y, np.clip(P, -1, 1), np.linspace(-1.0001, 1.0001, 21), cmap=cmap)
        ax.contour(X, Y, P, [0], colors="#1a1a19", linewidths=1.3)
        for a, b in d.E:
            ax.plot(*d.V[[a, b]].T, color="#52514e", lw=0.8, ls=(0, (2, 2)))
        tok = s[f"xyz_fps_{min(int(k[8:]) for k in s if k.startswith('xyz_fps_'))}"]
        ax.plot(*tok.T, ".", color="#8a8985", ms=1.2)
        ax.set_title(f"{d.name}   IoU {ious[i]:.3f}", fontsize=8, color="#1a1a19")
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#c9c8c4")
    for ax in axes.ravel()[len(idx):]:
        ax.axis("off")
    handles = [Line2D([], [], color="#1a1a19", lw=1.3, label="decoded zero level set"),
               Line2D([], [], color="#52514e", lw=0.8, ls=(0, (2, 2)), label="true boundary"),
               Line2D([], [], color="#8a8985", marker=".", ls="", label="tokens (xyz_fps_N)")]
    fig.suptitle(title, y=0.995, fontsize=9, color="#1a1a19")
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.465, 0.975), ncol=3,
               frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 0.93, 0.955))
    cax = fig.add_axes((0.94, 0.15, 0.012, 0.7))
    cb = fig.colorbar(cf, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label("decoded sdf / tau  (inside < 0)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.savefig(out, dpi=120, facecolor="#fcfcfb")
    plt.close(fig)


def save(path, model, cfg, step, metrics):
    """main_ae.py's layout: the model under `cond_vae.`, as the Geom-AE half of
    KLAutoEncoderJoint."""
    sd = {f"cond_vae.{k}": v.detach().cpu() for k, v in model.state_dict().items()}
    torch.save({"model": sd, "config": cfg, "step": step, "metrics": metrics}, path)


def load_geom_ae(path, device="cpu"):
    """Frozen Geom-AE for downstream use (SQ-VAE / SQ-Diffuse): returns (model, config)."""
    ck = torch.load(path, map_location=device, weights_only=False)
    model = build_model(ck["config"]).to(device)
    model.load_state_dict({k[len("cond_vae."):]: v for k, v in ck["model"].items() if k.startswith("cond_vae.")})
    model.requires_grad_(False)
    return model.eval(), ck["config"]


@torch.no_grad()
def encode_obj(model, cfg, path, device="cpu", seed=0):
    """Sec. 4 steps 1-3 for one OBJ: normalise, FPS N and 4N area samples, SDF features,
    then the geometry code from the posterior mean. Returns (code (1, N, dim), batch,
    Domain); `batch` holds xyz_fps_{N} (token positions), and Domain.centre/scale map
    them back to the OBJ's frame."""
    d = Domain(path)
    n = cfg["num_latents"]
    d.fps = fps_sets(d, (n, 4 * n), np.random.default_rng(seed))
    s = make_sample(d, dict(cfg, n_query=1, n_near=1, n_global=1), np.random.default_rng(seed), False)
    b = to_dev({k: v[None] for k, v in s.items() if "fps" in k}, device)
    return model.return_last_features_and_norm(model.encode(b, is_mode=True)["latent"]), b, d


def lr_at(step, cfg):
    """Linear warm-up, then half-cycle cosine to min_lr (the release's lr_sched, per step)."""
    warm = max(1, int(cfg["warmup_frac"] * cfg["steps"]))
    if step < warm:
        return cfg["lr"] * (step + 1) / warm
    t = (step - warm) / max(1, cfg["steps"] - warm)
    return cfg["min_lr"] + (cfg["lr"] - cfg["min_lr"]) * 0.5 * (1 + math.cos(math.pi * t))


def train(cfg, train_set, val_eval, train_eval, log):
    dev = pick_device()
    torch.manual_seed(cfg["seed"])
    model = build_model(cfg).to(dev)
    log(f"device={dev}  parameters={sum(p.numel() for p in model.parameters()):,}")
    opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    loader = torch.utils.data.DataLoader(
        DomainStream(train_set, cfg), batch_size=cfg["batch"], num_workers=cfg["workers"],
        persistent_workers=cfg["workers"] > 0, prefetch_factor=4 if cfg["workers"] > 0 else None)
    best, t0, run = float("inf"), time.time(), {}
    for step, batch in enumerate(loader, 1):
        for g in opt.param_groups:
            g["lr"] = lr_at(step - 1, cfg)
        loss, terms = calc_loss(model, to_dev(batch, dev), cfg["kl_weight"], is_mode=cfg["posterior"] == "mean")
        opt.zero_grad(set_to_none=True)
        loss.backward()
        gnorm = torch.nn.utils.clip_grad_norm_(model.parameters(), cfg["clip_grad"])
        opt.step()
        terms["loss"], terms["gnorm"] = loss, gnorm
        for k, v in terms.items():
            run[k] = run.get(k, 0.0) + v.detach().item()
        if step % cfg["log_every"] == 0 or step == 1:
            n = 1 if step == 1 else cfg["log_every"]
            el = time.time() - t0
            log(f"step {step:6d}  lr {opt.param_groups[0]['lr']:.2e}  "
                + "  ".join(f"{k} {v / n:.4f}" for k, v in run.items())
                + f"  ({el / step * 1000:.0f} ms/step, {el / 60:.1f} min)")
            run = {}
        if step % cfg["eval_every"] == 0 or step == cfg["steps"]:
            mv, mt = val_eval(model, dev), train_eval(model, dev)
            log(f"  eval {step:6d}  val band_l1 {mv['band_l1']:.4f} iou {mv['iou']:.4f} bnd {mv['bnd']:.4f}"
                f"  |  train band_l1 {mt['band_l1']:.4f} iou {mt['iou']:.4f} bnd {mt['bnd']:.4f}")
            save(os.path.join(cfg["out"], "last.pt"), model, cfg, step, {"val": mv, "train": mt})
            if mv["band_l1"] < best:
                best = mv["band_l1"]
                save(os.path.join(cfg["out"], "best.pt"), model, cfg, step, {"val": mv, "train": mt})
        if step >= cfg["steps"]:
            break
    return dev


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--data", nargs="+", default=DEFAULT_DIRS, help="directories of .obj domains")
    ap.add_argument("--out", default=os.path.join(REPO, "py", "squadgen", "runs", "geom_ae_2d"))
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--batch", type=int, default=8)                  # release: 16 per GPU; 8 fits 8 GB
    ap.add_argument("--n_tokens", type=int, default=512)             # N; keys/values 4N (Sec. 3.5)
    ap.add_argument("--n_query", type=int, default=1024)             # targets: area samples of Omega
    ap.add_argument("--n_near", type=int, default=1024)              # targets: near dOmega
    ap.add_argument("--n_global", type=int, default=512)             # targets: uniform in [-1, 1]^2
    ap.add_argument("--tau", type=float, default=0.1)                # SDF truncation (Sec. 2.7)
    ap.add_argument("--sigma_near", type=float, default=0.02)
    ap.add_argument("--w_query", type=float, default=1.0)            # Sec. 3.1 item 7: (1, 1, 0.1)
    ap.add_argument("--w_near", type=float, default=1.0)
    ap.add_argument("--w_global", type=float, default=0.1)           # 1.0 if it collapses (Sec. 3.2)
    ap.add_argument("--kl_weight", type=float, default=0.0)          # release: 0 for the Geom-AE
    ap.add_argument("--posterior", choices=("sample", "mean"), default="sample",
                    help="what the decoder sees in training (release: a posterior sample)")
    ap.add_argument("--dim", type=int, default=256)                  # Sec. 3.5: 256 / 8 / 32
    ap.add_argument("--depth", type=int, default=8)
    ap.add_argument("--latent_dim", type=int, default=32)
    ap.add_argument("--heads", type=int, default=8)                  # release hard-codes 8 x 64
    ap.add_argument("--dim_head", type=int, default=64)
    ap.add_argument("--drop_path_rate", type=float, default=0.1)     # release default
    ap.add_argument("--lr", type=float, default=2e-4)                # release: 1e-4
    ap.add_argument("--min_lr", type=float, default=1e-6)
    ap.add_argument("--warmup_frac", type=float, default=0.03)
    ap.add_argument("--wd", type=float, default=0.01)                # torch AdamW default, as the release
    ap.add_argument("--clip_grad", type=float, default=1.0)
    ap.add_argument("--no_aug", action="store_true")
    ap.add_argument("--val_frac", type=float, default=0.1)
    ap.add_argument("--workers", type=int, default=2)                # each holds torch + the corpus
    ap.add_argument("--log_every", type=int, default=100)
    ap.add_argument("--eval_every", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--eval_only", action="store_true", help="evaluate and plot best.pt")
    args = ap.parse_args()
    cfg = dict(vars(args), num_latents=args.n_tokens, augment=not args.no_aug)

    os.makedirs(args.out, exist_ok=True)
    logf = open(os.path.join(args.out, "log.txt"), "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  " + " ".join(f"{k}={v}" for k, v in sorted(cfg.items()) if k != "data"))
    # The token-set cache is a one-off CPU job before any model exists, so it may use more
    # processes than the DataLoader, which shares memory with the model.
    domains = load_domains(args.data, args.n_tokens, args.out, min(6, os.cpu_count() or 1))
    train_set, val_set = split_by_family(domains, args.val_frac, args.seed)
    log(f"{len(domains)} domains: {len(train_set)} train, {len(val_set)} val "
        f"({len({family(d.path) for d in val_set})} held-out families: "
        f"{' '.join(sorted({family(d.path) for d in val_set}))})")
    rng = np.random.default_rng(args.seed)
    val_eval = GridEval(val_set, cfg)
    train_eval = GridEval([train_set[i] for i in rng.choice(len(train_set), min(48, len(train_set)), replace=False)], cfg)

    dev = pick_device() if args.eval_only else train(cfg, train_set, val_eval, train_eval, log)
    model, _ = load_geom_ae(os.path.join(args.out, "best.pt"), dev)
    step = torch.load(os.path.join(args.out, "best.pt"), map_location="cpu", weights_only=False)["step"]
    for name, ev in (("val", val_eval), ("train", train_eval)):
        m, preds, ious = ev(model, dev, keep=True)
        log(f"best.pt (step {step}) on {name}: band_l1 {m['band_l1']:.4f}  iou {m['iou']:.4f}  bnd {m['bnd']:.4f}"
            f"  (IoU min {np.min(ious):.3f}, median {np.median(ious):.3f})")
        plot_reconstructions(ev, preds, ious, os.path.join(args.out, f"recon_{name}.png"),
                             f"Geom-AE 2D, best.pt (step {step}), {name} domains: worst and best IoU")
    log(f"plots in {args.out}")


if __name__ == "__main__":
    main()
