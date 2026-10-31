"""
geom_ae.py -- train SQuadGen's Geometry Autoencoder (Geom-AE) on planar domains.

This is the 2D Geom-AE of docs/squadgen_2d.md Sec. 3.2, trained on CrossGen domains (by default data/meshes/singlemat
and data/meshes/mambo: 442 single-material domains; any directories of OBJs with --data).

The network is a port of the release's KLAutoEncoder (SQuadGen squadgen/network/models_ae.py) with the Sec. 3.1 edits
for 2D: positions are 2-d (POS_DIM = 2) and the one point feature is the truncated SDF where the release has the
normal. Every module keeps the release's name and shape, and the state dict is saved with the `cond_vae.` prefix that
main_ae.py writes, so a checkpoint loads into the edited release class as the frozen geometry encoder of SQ-VAE
(train_sqvae.sh --pretrained_geom_vae_pth, which keeps the keys that start with `cond_vae.`). What is kept from the
release, and why:

* Encoder: N token points (xyz_fps_N) cross-attend to 4N key/value points (xyz_fps_4N) through one head as wide as
  the model, then a FeedForward, then a diagonal Gaussian bottleneck of `latent_dim` channels. The bottleneck is
  sampled in training with KL weight 0 (Sec. 1.3, "Geom-AE"), so it is a noisy narrow layer rather than a prior.
* proj back to `dim`, `depth` pre-norm self-attention blocks (8 heads x 64 whatever `dim` is, because
  create_vae_from_config hard-codes it; Sec. 3.5), and the decoder's context LayerNorm. Its output,
  `return_last_features_and_norm`, is the N x dim geometry code that SQ-VAE and SQ-Diffuse consume frozen. Token i
  sits at xyz_fps_N[i], the point SQ-VAE's query i reads (Sec. 3.2, "Token alignment").
* Decoder: one single-head cross-attention from the query point's Fourier embedding to the tokens, then a 3-layer
  MLP head per output key, zero-initialised.

What changes for 2D (Sec. 3.2, "2D", and the Appendix D loader):

* Tokens and keys/values are AREA samples of Omega, not boundary samples. SQ-VAE's query points are area samples too
  (the NPZ's fps sets lie on the split quad mesh, Sec. 2.4), and token i has to sit at query point i.
* The feature is sdf = clip(phi / tau, -1, 1), tau = 0.1, phi < 0 inside (Sec. 2.7 item 3). On area samples it lies
  in [-1, 0] and tells a token how far, up to tau, it is from dOmega. Domain shape is encoded both by where the
  points are and by these values.
* The target is the same truncated SDF at three groups of query points: output keys `sdf` (area samples, the
  release's "surface" group `xyz_query`), `sdf_near` (boundary samples + N(0, 0.02^2)) and `sdf_global` (uniform in
  [-1, 1]^2), weighted (1, 1, 0.1) (Sec. 3.1 item 7). The truncation is what keeps the 0.1 box weight from
  collapsing the output to phi = 0, which an untruncated SDF did (Sec. 3.2, "Collapse risk"); the release truncates
  its UDF at 0.08 for the same reason. --w_global 1 is the fallback.

Data (domains.py). Each domain's FPS token sets are computed once and cached in --out. Augmentation is the 2D one of
Sec. 2.7 item 1 (geometry.Deformation: rotation about z, reflection, anisotropic scale and a smooth warp), and a
deformation that folds a triangle is rejected. The inference normalisation of Sec. 4 step 1 follows (centre, then
0.999 / max radius), so training sees domains in the frame inference will hand it. The train/val split is by family
(the CAD model), and the held-out families are recorded in the checkpoint, where SQ-VAE and SQ-Diffuse take their own
split from.

Run
---
    python -m squadgen2d geom-ae                                          # train -> py/squadgen/runs/geom_ae_2d
    python -m squadgen2d geom-ae --out /tmp/smoke --steps 200             # smoke test
    python -m squadgen2d geom-ae --eval_only                              # plots from best.pt
    python -m squadgen2d geom-ae --data data/objs/2d_data --init py/squadgen/runs/geom_ae_2d/best.pt \\
        --out py/squadgen/runs/geom_ae_2d_more --steps 10000 --lr 1e-4    # more steps from a checkpoint

Continuing (--init). A checkpoint holds the weights but no optimizer state, so a continuation is a new run that
starts from them: a fresh AdamW with its own warm-up and cosine over --steps (step numbers start again at 1; the
config records init and init_step). The network's shape and tau come from the checkpoint, since the weights fit only
that network and predict SDF in units of that tau; every other flag is this run's. The loaded model is evaluated at
step 0 and saved as the first best.pt, so best.pt never ends worse than the start. When --out has no token-set cache,
the checkpoint directory's is copied (it is the first ten minutes to an hour of a large run). Keep --data, --val_frac
and --seed as they were: otherwise the held-out families change and validation counts domains the start was trained
on (logged as a note when they differ). The first run's cosine ended at min_lr, so a lower --lr (1e-4, 5e-5)
continues more gently than starting over at 2e-4, which first undoes some of what the end of that cosine refined.

Memory: on an 8 GB Apple M2, batch 8 with the default query counts takes ~3.3 GB of MPS memory at ~0.5 s/step; batch
12 needs 4.3 GB and batch 16 swapped (100 s/step). On an RTX 5090, 100k steps took 53 min.

Outputs in --out: fps_cache_N{N}.npz (token sets), best.pt / last.pt (checkpoint: model under `cond_vae.`, config,
step, metrics), recon_val.png / recon_train.png (decoded SDF against the true boundary), log.txt. Downstream code
calls load_geom_ae() and encode_obj().
"""

import argparse
import math
import os
import shutil
import time

import numpy as np
import torch
import torch.nn as nn

from . import geometry as geo
from . import training as tr
from .domains import (Domain, default_fps_workers, family, fps_sets, load_domains, obj_paths, split_by_family)
from .layers import KEY_DIM_DICT, Attention, FeedForward, PointEmbed, PreNorm, xavier_init
from .paths import CORPUS, DEFAULT_GEOM_AE, from_repo, portable

GEOM_KEYS = ["sdf"]


# =============================================================================
# 1. DATA
# =============================================================================


def make_sample(d, cfg, rng, aug):
    """One example under the release's batch keys (Appendix D, Geom-AE part only).
    Encoder input : xyz_fps_{N}, sdf_fps_{N} (tokens); xyz_fps_{4N}, sdf_fps_{4N} (keys/values).
    Targets       : sdf_query at xyz_query (area samples), sdf_near at xyz_near, and sdf_global at xyz_global, all
                    truncated SDF in units of tau."""
    V = geo.normalise(geo.augment(d.V, d.F, rng), d.F, d.bverts) if aug else d.V
    A, B = V[d.E[:, 0]], V[d.E[:, 1]]
    tau, out = cfg["tau"], {}
    for n in (cfg["n_tokens"], 4 * cfg["n_tokens"]):
        tri, bary = d.fps[n]
        c = rng.integers(len(tri))
        x = geo.point_from_bary(V, d.F, tri[c], bary[c].astype(np.float64))
        out[f"xyz_fps_{n}"] = x.astype(np.float32)
        out[f"sdf_fps_{n}"] = geo.truncated_sdf(x, A, B, tau, known_inside=True)
    x = geo.point_from_bary(V, d.F, *geo.sample_area(V, d.F, d.cdf, cfg["n_query"], rng))
    out["xyz_query"], out["sdf_query"] = x.astype(np.float32), geo.truncated_sdf(x, A, B, tau, True)
    L = np.linalg.norm(B - A, axis=1)
    ie = np.minimum(np.searchsorted(np.cumsum(L), rng.random(cfg["n_near"]) * L.sum()), len(L) - 1)
    x = (A[ie] + rng.random((len(ie), 1)) * (B[ie] - A[ie])
         + rng.normal(0, cfg["sigma_near"], (len(ie), 2)))
    out["xyz_near"], out["sdf_near"] = x.astype(np.float32), geo.truncated_sdf(x, A, B, tau)
    x = rng.uniform(-1, 1, (cfg["n_global"], 2))
    out["xyz_global"], out["sdf_global"] = x.astype(np.float32), geo.truncated_sdf(x, A, B, tau)
    return out


class DomainStream(torch.utils.data.IterableDataset):
    """Endless stream of augmented samples drawn uniformly over the training domains. Each DataLoader worker seeds
    its own generator from (seed, worker id)."""

    def __init__(self, domains, cfg):
        self.domains, self.cfg = domains, cfg

    def __iter__(self):
        info = torch.utils.data.get_worker_info()
        rng = np.random.default_rng([self.cfg["seed"], info.id if info is not None else 0])
        while True:
            d = self.domains[rng.integers(len(self.domains))]
            yield make_sample(d, self.cfg, rng, self.cfg["augment"])


# =============================================================================
# 2. NETWORK  (the release's KLAutoEncoder, models_ae.py, with the Sec. 3.1 edits)
# =============================================================================
# The release also loops the key/value points over "groups" and max-pools the encoder output across them. A [B, 4N, 2]
# input is one group, which makes that a no-op, so it is left out.


class KLAutoEncoder(nn.Module):
    """The release's KLAutoEncoder as the Geom-AE uses it: is_context_norm = 1 (so the encoder cross-attention
    normalises its context and the decoder normalises both sides), no learnable latents, no encoder self-attention,
    no decoder FeedForward, no QK-norm."""

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
        xavier_init(self)
        nn.init.normal_(self.point_embed.mlp_pos.weight, std=1)
        for k in self.in_keys_list_mlp:
            nn.init.normal_(getattr(self.point_embed, f"mlp_{k}")[0].weight, std=0.02)
        for k in self.out_keys_list_mlp:
            nn.init.constant_(self.to_outputs[k][-1].weight, 0)
            nn.init.constant_(self.to_outputs[k][-1].bias, 0)

    def encode(self, batch, is_mode=False):
        """-> {"latent": (B, N, latent_dim), "kl": (B,), "mean", "logvar"}. A sample of the posterior in training,
        its mean with is_mode (what SQ-VAE uses, use_geom_cond_mean=1)."""
        n = self.num_latents
        ctx = lambda s: {k: batch[f"{k}_fps_{s}"] for k in self.in_keys_list_mlp}
        x = self.point_embed(batch[f"xyz_fps_{n}"], ctx(n))
        pc = self.point_embed(batch[f"xyz_fps_{4 * n}"], ctx(4 * n))
        cross_attn, cross_ff = self.cross_attend_blocks
        x = cross_attn(x, context=pc) + x
        x = cross_ff(x) + x
        # DiagonalGaussianDistribution: logvar clamped to [-10, 10], KL averaged over tokens and channels (dims [1, 2]).
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
        """x from return_last_features (NOT the normed code: the decoder's PreNorm applies norm_context itself)."""
        h = self.decoder_cross_attn(self.point_embed.query(queries), context=x)
        return {k: self.to_outputs[k](h) for k in self.out_keys_list_mlp}

    def decode(self, z, queries):
        return self.query_last_features(self.return_last_features(z), queries)

    def decode_batch(self, z, batch):
        """Every output key at its own query group: sdf at xyz_query, sdf_near at xyz_near, sdf_global at
        xyz_global, decoded in one pass."""
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
    """The release's get_label: `<k>_near` / `<k>_global` are batch keys themselves, a surface key reads
    `<k>_query`."""
    return key if group_of(key) != "query" else f"{key}_query"


MODEL_KEYS = ("depth", "dim", "num_latents", "latent_dim", "heads", "dim_head", "drop_path_rate")


def build_model(cfg):
    return KLAutoEncoder(in_keys_list=GEOM_KEYS, out_keys_list=("sdf", "sdf_near", "sdf_global"),
                         out_weights_list=(cfg["w_query"], cfg["w_near"], cfg["w_global"]),
                         **{k: cfg[k] for k in MODEL_KEYS})


def save(path, model, cfg, step, metrics):
    """main_ae.py's layout: the model under `cond_vae.`, as the Geom-AE half of KLAutoEncoderJoint."""
    sd = {f"cond_vae.{k}": v.detach().cpu() for k, v in model.state_dict().items()}
    torch.save({"model": sd, "config": portable(cfg), "step": step, "metrics": metrics}, path)


def load_geom_ae(path, device="cpu"):
    """Frozen Geom-AE for downstream use (SQ-VAE / SQ-Diffuse): returns (model, config)."""
    ck = torch.load(path, map_location=device, weights_only=False)
    model = build_model(ck["config"]).to(device)
    model.load_state_dict({k[len("cond_vae."):]: v for k, v in ck["model"].items() if k.startswith("cond_vae.")})
    model.requires_grad_(False)
    return model.eval(), ck["config"]


@torch.no_grad()
def encode_obj(model, cfg, path, device="cpu", seed=0):
    """Sec. 4 steps 1-3 for one OBJ: normalise, FPS N and 4N area samples, SDF features, then the geometry code from
    the posterior mean. Returns (code (1, N, dim), batch, Domain); `batch` holds xyz_fps_{N} (token positions), and
    Domain.centre/scale map them back to the OBJ's frame."""
    d = Domain(path)
    n = cfg["num_latents"]
    d.fps = fps_sets(d, (n, 4 * n), np.random.default_rng(seed))
    s = make_sample(d, dict(cfg, n_query=1, n_near=1, n_global=1), np.random.default_rng(seed), False)
    b = tr.to_dev({k: v[None] for k, v in s.items() if "fps" in k}, device)
    return model.return_last_features_and_norm(model.encode(b, is_mode=True)["latent"]), b, d


# =============================================================================
# 3. EVALUATION
# =============================================================================


def calc_loss(model, b, kl_weight, is_mode=False):
    """engine_ae.calc_loss for the Geom-AE alone: mean L1 per output key (per sample over points and channels, then
    over the batch), weighted, plus kl_weight x KL."""
    enc = model.encode(b, is_mode=is_mode)
    pred = model.decode_batch(enc["latent"], b)
    terms = {k: (pred[k] - b[label_of(k)]).abs().mean(dim=[1, 2]).mean() for k in model.out_keys_list}
    terms["kl"] = enc["kl"].mean()
    loss = sum(w * terms[k] for k, w in zip(model.out_keys_list, model.out_weights_list))
    return loss + kl_weight * terms["kl"], terms


class GridEval:
    """Fixed, unaugmented evaluation of a domain set, using the posterior mean. The decoded SDF on a res^2 grid of
    [-1, 1]^2 is compared against the exact one, and three numbers are reported:
      band_l1  L1 error where the true truncated SDF is inside the band (|phi| < tau), in units of tau. That is
               where the target carries information.
      iou      IoU of the decoded interior {phi_hat < 0} with Omega.
      bnd      mean |phi_hat| on dOmega in units of tau: roughly how far the decoded zero level set sits from the
               true boundary, as a fraction of tau."""

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
            self.items.append((d, s, geo.truncated_sdf(self.grid, A, B, cfg["tau"])[:, 0], on.astype(np.float32)))
        self.grid = self.grid.astype(np.float32)

    @torch.no_grad()
    def __call__(self, model, dev, keep=False):
        was = model.training
        model.eval()
        band, iou, bnd, preds = [], [], [], []
        for d, s, tgt, on in self.items:
            b = tr.to_dev({k: v[None] for k, v in s.items()}, dev)
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
    """A panel per domain, worst IoU first, so that the failures are the first thing seen. Each panel shows the
    decoded truncated SDF as filled colour, its zero level set as a solid line, the true dOmega dashed, and the
    tokens as grey dots."""
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
        # The decoder is not clipped and overshoots +-1 a little where the target saturates; contourf would leave
        # those regions unfilled, which reads as an artefact.
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


# =============================================================================
# 4. TRAINING
# =============================================================================


def train(cfg, train_set, val_eval, train_eval, log, init=None):
    """init: a checkpoint dict to continue from (see Continuing in the module docstring)."""
    dev = cfg["device"]
    amp = tr.amp_dtype(dev, cfg["amp"])
    torch.manual_seed(cfg["seed"])
    model = build_model(cfg).to(dev)
    log(f"device={dev}  amp={amp}  workers={cfg['workers']}  parameters={sum(p.numel() for p in model.parameters()):,}")
    best = float("inf")
    if init is not None:
        model.load_state_dict({k[len("cond_vae."):]: v for k, v in init["model"].items() if k.startswith("cond_vae.")})
        # The start is the best so far, so a continuation that only gets worse leaves best.pt as it began.
        mv, mt = val_eval(model, dev), train_eval(model, dev)
        log(f"  eval {0:6d}  val band_l1 {mv['band_l1']:.4f} iou {mv['iou']:.4f} bnd {mv['bnd']:.4f}"
            f"  |  train band_l1 {mt['band_l1']:.4f} iou {mt['iou']:.4f} bnd {mt['bnd']:.4f}"
            f"  (start: {cfg['init']}, step {init['step']} of its run)")
        best = mv["band_l1"]
        save(os.path.join(cfg["out"], "best.pt"), model, cfg, 0, {"val": mv, "train": mt})
        tr.empty_cache(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    scaler = tr.grad_scaler(dev, amp)
    loader = tr.make_loader(DomainStream(train_set, cfg), cfg)
    t0, run = time.time(), {}
    for step, batch in enumerate(loader, 1):
        for g in opt.param_groups:
            g["lr"] = tr.lr_at(step - 1, cfg)
        with tr.autocast(dev, amp):
            loss, terms = calc_loss(model, tr.to_dev(batch, dev), cfg["kl_weight"], is_mode=cfg["posterior"] == "mean")
        gnorm = tr.train_step(model, opt, scaler, loss, cfg["clip_grad"])
        terms["loss"], terms["gnorm"] = loss, gnorm
        tr.accumulate(run, terms)
        if step % cfg["log_every"] == 0 or step == 1:
            el = time.time() - t0
            log(f"step {step:6d}  lr {opt.param_groups[0]['lr']:.2e}  {tr.averages(run)}"
                + f"  ({el / step * 1000:.0f} ms/step, {el / 60:.1f} min{tr.memory_note(dev)})")
            run = {}
        if step % cfg["eval_every"] == 0 or step == cfg["steps"]:
            mv, mt = val_eval(model, dev), train_eval(model, dev)
            log(f"  eval {step:6d}  val band_l1 {mv['band_l1']:.4f} iou {mv['iou']:.4f} bnd {mv['bnd']:.4f}"
                f"  |  train band_l1 {mt['band_l1']:.4f} iou {mt['iou']:.4f} bnd {mt['bnd']:.4f}")
            save(os.path.join(cfg["out"], "last.pt"), model, cfg, step, {"val": mv, "train": mt})
            if mv["band_l1"] < best:
                best = mv["band_l1"]
                save(os.path.join(cfg["out"], "best.pt"), model, cfg, step, {"val": mv, "train": mt})
            tr.empty_cache(dev)
        if step >= cfg["steps"]:
            break


def main(argv=None):
    ap = argparse.ArgumentParser(prog="squadgen2d geom-ae", description=__doc__.split("\n")[1].strip())
    ap.add_argument("--data", nargs="+", default=CORPUS, help="directories of .obj domains, searched recursively")
    ap.add_argument("--out", default=DEFAULT_GEOM_AE)
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
    ap.add_argument("--fps_workers", type=int,
                    help="processes for the one-off token-set cache (default: every core on CUDA, 6 elsewhere)")
    tr.add_device_args(ap)
    ap.add_argument("--eval_only", action="store_true", help="evaluate and plot best.pt")
    ap.add_argument("--init", help="continue from this Geom-AE checkpoint (best.pt / last.pt; see Continuing)")
    args = ap.parse_args(argv)
    init = None
    if args.init:
        args.init = args.init if os.path.exists(args.init) else from_repo(args.init)
        init = torch.load(args.init, map_location="cpu", weights_only=False)
        for k in ("depth", "dim", "latent_dim", "heads", "dim_head", "drop_path_rate", "n_tokens", "tau"):
            setattr(args, k, init["config"][k])   # the network the weights fit, and the units they predict in
    dev = tr.pick_device(args.device)
    cfg = dict(vars(args), num_latents=args.n_tokens, augment=not args.no_aug, device=dev,
               workers=tr.default_workers(dev) if args.workers is None else args.workers,
               init_step=init["step"] if init else None)

    log = tr.RunLog(args.out)
    log.header(cfg, skip=("data",))
    if init:
        name = f"fps_cache_N{args.n_tokens}.npz"
        src, dst = os.path.join(os.path.dirname(args.init), name), os.path.join(args.out, name)
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.copyfile(src, dst)
            log(f"token-set cache copied from {os.path.dirname(args.init)}")
    # The token-set cache is a one-off CPU job before any model exists, so it may use more processes than the
    # DataLoader, which shares memory with the model.
    domains = load_domains(args.data, args.n_tokens, args.out, args.fps_workers or default_fps_workers(dev))
    train_set, val_set = split_by_family(domains, args.val_frac, args.seed)
    # Recorded in every checkpoint: SQ-VAE and SQ-Diffuse validate on these families (domains.geom_ae_val_families).
    cfg["val_families"] = sorted({family(d.path) for d in val_set})
    log(f"{len(domains)} domains: {len(train_set)} train, {len(val_set)} val "
        f"({len(cfg['val_families'])} held-out families: {' '.join(cfg['val_families'])})")
    if init:
        # The split is a function of the domain names, --val_frac and --seed; compared by name, so the same domains
        # passed as one parent directory or as its subdirectories count as the same.
        ic = init["config"]
        before = {os.path.basename(p) for p in obj_paths([from_repo(d) for d in ic["data"]])}
        if (before != {os.path.basename(d.path) for d in domains}
                or (ic["val_frac"], ic["seed"]) != (args.val_frac, args.seed)):
            log("note: the domains, --val_frac or --seed differ from the checkpoint's run (or its --data is not "
                "here), so the held-out families may differ: validation may count domains the start trained on")
    rng = np.random.default_rng(args.seed)
    val_eval = GridEval(val_set, cfg)
    train_eval = GridEval([train_set[i] for i in rng.choice(len(train_set), min(48, len(train_set)), replace=False)], cfg)

    if not args.eval_only:
        train(cfg, train_set, val_eval, train_eval, log, init)
    model, _ = load_geom_ae(os.path.join(args.out, "best.pt"), dev)
    step = torch.load(os.path.join(args.out, "best.pt"), map_location="cpu", weights_only=False)["step"]
    for name, ev in (("val", val_eval), ("train", train_eval)):
        m, preds, ious = ev(model, dev, keep=True)
        log(f"best.pt (step {step}) on {name}: band_l1 {m['band_l1']:.4f}  iou {m['iou']:.4f}  bnd {m['bnd']:.4f}"
            f"  (IoU min {np.min(ious):.3f}, median {np.median(ious):.3f})")
        plot_reconstructions(ev, preds, ious, os.path.join(args.out, f"recon_{name}.png"),
                             f"Geom-AE 2D, best.pt (step {step}), {name} domains: worst and best IoU")
    log(f"plots in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
