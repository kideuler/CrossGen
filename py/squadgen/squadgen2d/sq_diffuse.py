"""
sq_diffuse.py -- train SQuadGen's SQ-Diffuse on CrossGen layouts of planar domains.

This is the SQ-Diffuse of docs/squadgen_2d.md Sec. 3.4. It learns to generate SQ-VAE latents from the geometry alone:
given the frozen Geom-AE's code of a domain, an N x latent_dim latent whose decoded fields QuadExtraction turns into a
layout of that domain (Sec. 1.2, the last line of the data flow). The SQ-VAE is sq_vae.py's best.pt and the Geom-AE
is the one inside it, both frozen.

The network is the release's SiT (SQuadGen squadgen/network/models_sit.py) as sqdiffuse.yaml configures it, with the
Sec. 3.1 edits for 2D (items 3-5: POS_DIM in PointEmbed, sdf in place of the normal). Module names and shapes are the
release's, so the state dict loads into the edited release class. What is kept, and why:

* Tokens. Token i is SQ-VAE latent token i, which sits at xyz_fps_N[i], the point Geom-AE token i was computed at
  (Sec. 3.2, "Token alignment"). The tokens are a set: SiT adds no positional embedding of its own, and the geometry
  is the only thing that tells token i where it is.
* Geometry conditioning before every block (is_last_feature = is_pecode_latent = is_pecode_latent_pe = 1): x +=
  Linear_l(code_i) + PointEmbed_l(xyz_i, sdf_i), with code the Geom-AE's return_last_features_and_norm output
  (in_channels_cond = its width) and a separate Linear and PointEmbed per block (Sec. 1.3, "Geometry conditioning").
  The paper describes a single addition at the input; the release does this.
* Time enters only through adaLN-Zero modulation, from a sinusoidal embedding of t in [0, 1]. Attention has QK-norm by
  nn.RMSNorm (qk_norm = 1); no learn_sigma. The adaLN layers and the final layer are zero-initialised, so the network
  starts as v = 0.
* No classifier-free guidance. The release's drop_prob only drops condition_global, which sqdiffuse.yaml leaves
  empty, so it does nothing (Sec. 1.3, "drop_prob").

Transport (Sec. 3.4, the release's transport/ with path ICPlan, prediction velocity, loss weight None): x_t = t x1 +
(1 - t) x0 with x0 ~ N(0, I) at t = 0 and the data at t = 1, t ~ U(0, 1) per sample, and the network regresses the
velocity x1 - x0 by MSE, averaged over tokens and channels. The data x1 are SQ-VAE posterior MEANS (is_mode = 1) of
augmented layouts, encoded on the fly in the main process exactly as the release does, normalised per channel by
target_mu / target_std from sq_vae.py's latent_stats.json. The statistics have to be the 2D SQ-VAE's own (Sec. 7,
"Recompute target_mu/target_std"); they are written into the config as the release's sqdiffuse.yaml holds them, and
checked against the SQ-VAE checkpoint's step.

Sampling is the release's deterministic ODE from t = 0 to 1: dopri5 with atol 1e-6, rtol 1e-3
(util.sample_with_spatial_smoothing -> Sampler.sample_ode). torchdiffeq is not a dependency here, so dopri5() below is
the Dormand-Prince 5(4) pair with torchdiffeq's initial step and step controller (one step size for the whole batch,
RMS error norm over every element). Inference smooths the initial noise first: Taubin smoothing over each token's 32
nearest tokens (Appendix C), with pairs whose segment crosses dOmega blocked, which in 2D does the job the release's
0.1-scaled normals do on a surface. The release trains on white noise and samples from smoothed noise (std ~0.56), so
the final evaluation reports both.

Evaluation. A domain has up to five layouts in the data (one per method), and SQ-Diffuse is meant to produce some good
layout, not a particular one. So a sample is scored against the CLOSEST of its domain's layouts: wall IoU, as
sq_vae.py defines it (cdf <= 0.1 decoded vs true), at 4096 area samples of Omega. The domain's geometry comes from its
OBJ exactly as inference will see it (geom_ae.encode_obj: Sec. 4 steps 1-3), and each layout's true fields are read
off its split mesh at the same points. The SQ-VAE's own reconstruction of every layout, scored the same way, is the
ceiling: SQ-Diffuse produces latents, and the decoder can do no better than it does on the encoder's. Model selection
uses the flow loss on held-out domains at fixed (t, x0) draws, which is deterministic; the sample IoU is logged beside
it with a fixed-step Heun sampler, which keeps the training-time numbers comparable between evaluations.

Training follows the release (train_sqdiffuse.sh): AdamW (weight decay 0.01, torch's default), lr 1e-4, gradient
clip 1, and on an NVIDIA GPU bf16 autocast (--amp, training.amp_dtype; fp32 on MPS, where it was slower). Only the SiT
runs under autocast. The frozen encoders that make its targets stay in fp32, unlike the release: bf16 rounds a latent
by ~0.4 %, about the SQ-VAE's own posterior spread, and the decoder cannot absorb errors that size. The release holds
the rate constant (--const_lr); the default here warms up and decays it by a half cosine, because the SQ-VAE run
gained most of its wall IoU (0.70 -> 0.78) only as its rate decayed. The default size is Sec. 3.5's 12 blocks x 6
heads x 64 (hidden 384), with the 512 tokens the Geom-AE and SQ-VAE were trained at.

Portability: the checkpoint records the SQ-VAE it was trained on by its path relative to the repository and by a hash
of its weights (vae_fingerprint), and holds the SiT's weights on the CPU, so a model trained on a cloud GPU loads on a
Mac once both run directories are copied to the same places under py/squadgen/runs; load_sq_diffuse() refuses an
SQ-VAE whose weights differ.

Memory: on an 8 GB Apple M2, batch 8 takes 4.6 GiB of MPS memory at ~1.15 s/step (SiT 0.95 s, the frozen Geom-AE and
SQ-VAE encoder 0.2 s), so 20k steps take ~6.5 h plus ~1 min per evaluation. Batch 16 needs 7.5 GiB on its own and
swapped (~160 s/step); in the full script it went NaN within 20 steps under that pressure, though the same steps alone
stayed finite.

Data: the field pack of a cdf-data directory (cdf_pack.py), as sq_vae.py reads it, with the SQ-VAE's split (the
families its Geom-AE held out, recorded in the SQ-VAE checkpoint). Samples are drawn uniformly over the training
layouts and augmented as there. The evaluation also reads the layouts' NPZs and their domains' OBJs. Sec. 6 step 6
asks first for an overfit: --overfit B13_1 trains and evaluates on that domain's layouts alone.

Run
---
    python -m squadgen2d sq-diffuse                                          # train -> py/squadgen/runs/sq_diffuse_2d
    python -m squadgen2d sq-diffuse --overfit B13_1 --steps 3000 \\
        --out py/squadgen/runs/sq_diffuse_2d_overfit                         # Sec. 6 step 6: overfit one domain
    python -m squadgen2d sq-diffuse --eval_only                              # metrics and plots of best.pt

Outputs in --out:

* best.pt / last.pt: the SiT state dict, its config (with target_mu / target_std and the SQ-VAE checkpoint it belongs
  to), step and metrics;
* samples_val.png / samples_train.png: per domain, the closest true layout's CDF and n_gen decoded samples;
* samples.csv: every final sample's wall IoU and closest layout;
* log.txt.

Downstream code (Sec. 4 inference) calls load_sq_diffuse(), domain_condition() and sample_latents(), then
sq_vae.decode_fields() and roundtrip's extraction functions.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import time

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from . import geometry as geo
from . import training as tr
from .cdf_pack import make_sample, open_pack
from .domains import Domain, geom_ae_val_families
from .geom_ae import GEOM_KEYS
from .geom_ae import make_sample as geom_sample
from .layers import PointEmbed
from .layout import METHODS, Layout, chart_sides, layout_fields
from .paths import DEFAULT_CDF, DEFAULT_SQ_DIFFUSE, DEFAULT_SQ_VAE, REPO, from_repo, portable, to_repo
from .roundtrip import extraction_mesh, locate
from .sq_vae import FIELD_RAMP, LayoutStream, decode_fields, geom_code, load_sq_vae, split_layouts

DEFAULT_VAE = os.path.join(DEFAULT_SQ_VAE, "best.pt")
EVAL_POINTS = 4096
# Taubin's lambda / mu as util.spatial_smoothing has them (Sec. 1.3, "Taubin weights").
TAUBIN = (0.4507499669, -0.4720265626)


# =============================================================================
# 1. NETWORK  (the release's SiT, models_sit.py, with the Sec. 3.1 edits)
# =============================================================================


def modulate(x, shift, scale):
    return x * (1 + scale.unsqueeze(1)) + shift.unsqueeze(1)


class TimestepEmbedder(nn.Module):
    """Sinusoidal features of t (t in [0, 1] as it is, no rescaling to 1000 steps), then a
    two-layer MLP."""

    def __init__(self, hidden_size, frequency_embedding_size=256, max_period=10000):
        super().__init__()
        self.mlp = nn.Sequential(nn.Linear(frequency_embedding_size, hidden_size), nn.SiLU(),
                                 nn.Linear(hidden_size, hidden_size))
        self.frequency_embedding_size = frequency_embedding_size
        self.max_period = max_period

    def forward(self, t):
        half = self.frequency_embedding_size // 2
        freqs = torch.exp(-math.log(self.max_period) * torch.arange(half, dtype=torch.float32, device=t.device) / half)
        args = t[:, None].float() * freqs[None]
        return self.mlp(torch.cat([torch.cos(args), torch.sin(args)], dim=-1))


class Attention(nn.Module):
    """timm's Attention as the release subclasses it: fused qkv with bias, per-head QK-norm
    (nn.RMSNorm over the head dimension), scaled-dot-product attention, output projection."""

    def __init__(self, dim, num_heads, qk_norm=True):
        super().__init__()
        self.num_heads, self.head_dim = num_heads, dim // num_heads
        self.qkv = nn.Linear(dim, dim * 3, bias=True)
        self.q_norm = nn.RMSNorm(self.head_dim) if qk_norm else nn.Identity()
        self.k_norm = nn.RMSNorm(self.head_dim) if qk_norm else nn.Identity()
        self.proj = nn.Linear(dim, dim)

    def forward(self, x):
        B, N, C = x.shape
        q, k, v = self.qkv(x).reshape(B, N, 3, self.num_heads, self.head_dim).permute(2, 0, 3, 1, 4).unbind(0)
        x = F.scaled_dot_product_attention(self.q_norm(q), self.k_norm(k), v)
        return self.proj(x.transpose(1, 2).reshape(B, N, C))


class Mlp(nn.Module):
    """timm's Mlp with the tanh GELU, under its module names (fc1, fc2)."""

    def __init__(self, dim, hidden):
        super().__init__()
        self.fc1, self.act, self.fc2 = nn.Linear(dim, hidden), nn.GELU(approximate="tanh"), nn.Linear(hidden, dim)

    def forward(self, x):
        return self.fc2(self.act(self.fc1(x)))


class SiTBlock(nn.Module):
    """Pre-norm attention and MLP, each shifted, scaled and gated by the time embedding
    (adaLN-Zero)."""

    def __init__(self, hidden_size, num_heads, mlp_ratio=4.0, qk_norm=True):
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_size, elementwise_affine=False, eps=1e-6)
        self.attn = Attention(hidden_size, num_heads, qk_norm)
        self.norm2 = nn.LayerNorm(hidden_size, elementwise_affine=False, eps=1e-6)
        self.mlp = Mlp(hidden_size, int(hidden_size * mlp_ratio))
        self.adaLN_modulation = nn.Sequential(nn.SiLU(), nn.Linear(hidden_size, 6 * hidden_size))

    def forward(self, x, c):
        shift_msa, scale_msa, gate_msa, shift_mlp, scale_mlp, gate_mlp = self.adaLN_modulation(c).chunk(6, dim=1)
        x = x + gate_msa.unsqueeze(1) * self.attn(modulate(self.norm1(x), shift_msa, scale_msa))
        return x + gate_mlp.unsqueeze(1) * self.mlp(modulate(self.norm2(x), shift_mlp, scale_mlp))


class FinalLayer(nn.Module):
    def __init__(self, hidden_size, out_channels):
        super().__init__()
        self.norm_final = nn.LayerNorm(hidden_size, elementwise_affine=False, eps=1e-6)
        self.linear = nn.Linear(hidden_size, out_channels)
        self.adaLN_modulation = nn.Sequential(nn.SiLU(), nn.Linear(hidden_size, 2 * hidden_size))

    def forward(self, x, c):
        shift, scale = self.adaLN_modulation(c).chunk(2, dim=1)
        return self.linear(modulate(self.norm_final(x), shift, scale))


class SiT(nn.Module):
    """The release's SiT with sqdiffuse.yaml's switches: is_last_feature = 1, qk_norm = 1,
    is_pecode_latent = 1, is_pecode_latent_pe = 1, condition_global_keys_list = [], no
    learn_sigma. The latent statistics are buffers that the state dict leaves out
    (persistent=False), as the release keeps them in the config, not the weights.

    `cond` is {"code": (B, N, in_channels_cond), "xyz": (B, N, 2), "sdf": (B, N, 1)}: the
    release's condition_params["pe_latent"] as {"code", "raw_pe": {"input", "context_dict"}}."""

    def __init__(self, in_channels=16, depth=12, n_heads=6, d_head=64, mlp_ratio=4.0, in_channels_cond=256,
                 qk_norm=1, target_mu=None, target_std=None):
        super().__init__()
        hidden = n_heads * d_head
        self.in_channels = self.channels = in_channels
        self.t_embedder = TimestepEmbedder(hidden)
        self.input_proj = nn.Linear(in_channels, hidden)
        self.blocks = nn.ModuleList([SiTBlock(hidden, n_heads, mlp_ratio, bool(qk_norm)) for _ in range(depth)])
        self.final_layer = FinalLayer(hidden, in_channels)
        self.pe_embed_list = nn.ModuleList([nn.Linear(in_channels_cond, hidden, bias=False) for _ in range(depth)])
        self.pe_embed_list2 = nn.ModuleList([PointEmbed(hidden, GEOM_KEYS) for _ in range(depth)])
        mu = torch.zeros(in_channels) if target_mu is None else torch.tensor(target_mu, dtype=torch.float32)
        std = torch.ones(in_channels) if target_std is None else torch.tensor(target_std, dtype=torch.float32)
        assert mu.shape == std.shape == (in_channels,), "target_mu / target_std need one entry per latent channel"
        self.register_buffer("target_mu", mu, persistent=False)
        self.register_buffer("target_std", std, persistent=False)
        self.initialize_weights()

    def initialize_weights(self):
        def basic(m):
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
        self.apply(basic)
        nn.init.normal_(self.t_embedder.mlp[0].weight, std=0.02)
        nn.init.normal_(self.t_embedder.mlp[2].weight, std=0.02)
        for block in self.blocks:
            nn.init.constant_(block.adaLN_modulation[-1].weight, 0)
            nn.init.constant_(block.adaLN_modulation[-1].bias, 0)
        nn.init.constant_(self.final_layer.adaLN_modulation[-1].weight, 0)
        nn.init.constant_(self.final_layer.adaLN_modulation[-1].bias, 0)
        nn.init.constant_(self.final_layer.linear.weight, 0)
        nn.init.constant_(self.final_layer.linear.bias, 0)

    def normalize(self, z):
        return (z - self.target_mu) / self.target_std

    def denormalize(self, z):
        return z * self.target_std + self.target_mu

    def cond_embed(self, cond):
        """What each block adds to the tokens before it. It depends on the geometry only, so a
        sampler computes it once per trajectory instead of once per function evaluation."""
        ctx = {k: cond[k] for k in GEOM_KEYS}
        return [lin(cond["code"]) + pe(cond["xyz"], ctx) for lin, pe in zip(self.pe_embed_list, self.pe_embed_list2)]

    def forward(self, x, t, cond, pe=None):
        """The velocity at normalised latents x (B, N, C) and times t (B,)."""
        pe = self.cond_embed(cond) if pe is None else pe
        x = self.input_proj(x)
        c = self.t_embedder(t)
        for block, p in zip(self.blocks, pe):
            x = block(x + p, c)
        return self.final_layer(x, c)


MODEL_KEYS = ("in_channels", "depth", "n_heads", "d_head", "mlp_ratio", "in_channels_cond", "qk_norm",
              "target_mu", "target_std")


def build_model(cfg):
    return SiT(**{k: cfg[k] for k in MODEL_KEYS})


def flow_loss(sit, x1, cond, t=None, x0=None):
    """transport.training_losses for the linear path and velocity prediction, per sample:
    mean over tokens and channels of (v(x_t, t) - (x1 - x0))^2. x1 is normalised. The release
    draws x0 first and then t, which is kept so that a seeded CPU run draws the same numbers."""
    x0 = torch.randn_like(x1) if x0 is None else x0
    t = torch.rand(x1.shape[0]).to(x1) if t is None else t
    tt = t[:, None, None]
    v = sit(tt * x1 + (1 - tt) * x0, t, cond)
    return ((v - (x1 - x0)) ** 2).mean(dim=[1, 2])


@torch.no_grad()
def encode_batch(vae, geom, b):
    """The release's get_latent_and_condition with is_mode = 1: the SQ-VAE posterior mean of a
    batch of layouts (B, N, latent_dim), not normalised, and the SiT condition at its tokens."""
    code = geom_code(geom, b)
    z = vae.encode(b, code, is_mode=True)["latent"]
    n = vae.num_latents
    return z, {"code": code, "xyz": b[f"xyz_fps_{n}"], "sdf": b[f"sdf_fps_{n}"]}


# =============================================================================
# 2. SAMPLING
# =============================================================================


def segments_cross(p, q, a, b, eps=1e-12):
    """Whether segment p-q properly crosses any of the segments a-b: p, q (..., 2), a, b (S, 2)
    -> (...) bool (Appendix C)."""
    def orient(u, v, w):
        return (v[..., 0] - u[..., 0]) * (w[..., 1] - u[..., 1]) - (v[..., 1] - u[..., 1]) * (w[..., 0] - u[..., 0])
    p, q = p[..., None, :], q[..., None, :]
    return (((orient(p, q, a) * orient(p, q, b)) < -eps) & ((orient(a, b, p) * orient(a, b, q)) < -eps)).any(-1)


def taubin_smooth(z, pts, a, b, k=32, ring_size=2, iters=10):
    """Appendix C for one domain: Taubin smoothing of noise z (N, C) over the k nearest tokens
    (the token itself included) with Gaussian weights whose width is ring_size x the distance to
    the nearest other token, as util.spatial_smoothing has it. In place of the release's normal
    term, a pair whose segment crosses dOmega (segments a-b) gets weight 0, so noise is not
    correlated across a slot or a thin wall. Ten alternating lambda / mu passes."""
    N = len(pts)
    k = min(k, N)
    sq = torch.cdist(pts, pts) ** 2
    sq_k, idx = torch.topk(sq, k=k, dim=1, largest=False, sorted=True)
    w = torch.exp(-sq_k / (2 * ring_size ** 2 * sq_k[:, 1:2].clamp_min(1e-20)))
    step = max(1, (1 << 22) // max(1, k * len(a)))
    blocked = torch.cat([segments_cross(pts[s:s + step, None].expand(-1, k, 2), pts[idx[s:s + step]], a, b)
                         for s in range(0, N, step)])
    w = torch.where(blocked, torch.zeros_like(w), w)
    w = w / w.sum(1, keepdim=True)
    for i in range(iters):
        c = TAUBIN[i % 2]
        z = c * (z[idx] * w[..., None]).sum(1) + (1 - c) * z
    return z


# Dormand-Prince 5(4): nodes, stages (the last row is the 5th-order solution, so the seventh
# stage is f at the new point and is reused as the next step's first, FSAL), and the
# difference between the 5th- and 4th-order weights for the error estimate.
DOPRI5_C = (1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0)
DOPRI5_A = ((1 / 5,), (3 / 40, 9 / 40), (44 / 45, -56 / 15, 32 / 9),
            (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
            (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
            (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84))
DOPRI5_E = (35 / 384 - 5179 / 57600, 0.0, 500 / 1113 - 7571 / 16695, 125 / 192 - 393 / 640,
            -2187 / 6784 + 92097 / 339200, 11 / 84 - 187 / 2100, -1 / 40)


def dopri5(f, x, t0=0.0, t1=1.0, atol=1e-6, rtol=1e-3, max_steps=10000):
    """Adaptive Dormand-Prince from t0 to t1 as torchdiffeq's 'dopri5' takes it: its initial
    step (_select_initial_step, order 4), error ratio (RMS over every element of
    err / (atol + rtol max(|x|, |x_new|))) and controller (safety 0.9, growth at most 10, shrink at
    least 0.2, exponent 1/5). The last step is cut to land on t1, where torchdiffeq would step
    past and interpolate. Returns (x(t1), function evaluations)."""
    def rms(v):
        return v.pow(2).mean().sqrt().item()

    k1 = f(x, t0)
    scale = atol + x.abs() * rtol
    d0, d1 = rms(x / scale), rms(k1 / scale)
    h0 = 1e-6 if d0 < 1e-5 or d1 < 1e-5 else 0.01 * d0 / d1
    d2 = rms((f(x + h0 * k1, t0 + h0) - k1) / scale) / h0
    h = min(100 * h0, max(1e-6, h0 * 1e-3) if max(d1, d2) <= 1e-15 else (0.01 / max(d1, d2)) ** 0.2)
    t, nfe = t0, 2
    for _ in range(max_steps):
        if t1 - t <= 1e-12:
            return x, nfe
        h = min(h, t1 - t)
        ks = [k1]
        for c, a in zip(DOPRI5_C, DOPRI5_A):
            y = x + h * sum(ai * ki for ai, ki in zip(a, ks) if ai != 0)
            ks.append(f(y, t + c * h))
        nfe += 6
        err = h * sum(e * k for e, k in zip(DOPRI5_E, ks) if e != 0)
        ratio = rms(err / (atol + rtol * torch.maximum(x.abs(), y.abs())))
        if ratio <= 1:
            t, x, k1 = t + h, y, ks[-1]
        h *= 10.0 if ratio == 0 else min(10.0, max(0.9 * ratio ** -0.2, 0.2 if ratio > 1 else 1.0))
    raise RuntimeError(f"dopri5: no convergence in {max_steps} steps (t = {t:.4f}, h = {h:.2e})")


def fixed_step(f, x, steps, method="heun"):
    """Euler or Heun (trapezoidal predictor-corrector) on a uniform grid of `steps` steps."""
    h = 1.0 / steps
    for i in range(steps):
        v = f(x, i * h)
        if method == "heun":
            v = 0.5 * (v + f(x + h * v, (i + 1) * h))
        x = x + h * v
    return x, steps * (2 if method == "heun" else 1)


@torch.no_grad()
def sample_latents(sit, cond, seeds, noise="taubin", bnd=None, method="dopri5", steps=12, atol=1e-6, rtol=1e-3):
    """Latents (B, N, C), denormalised, for B conditions at once: one initial noise per seed
    (a CPU generator each, so a seed gives the same sample on any device, as the release's
    StackedRandomGenerator does per device), Taubin-smoothed when noise == 'taubin' with each
    domain's boundary segments bnd[i] = (a, b) in the token frame, then the ODE from t = 0 to 1.
    Returns (latents, function evaluations)."""
    B, N = cond["xyz"].shape[:2]
    dev = cond["xyz"].device
    z = torch.stack([torch.randn((N, sit.channels), generator=torch.Generator().manual_seed(int(s))) for s in seeds])
    if noise == "taubin":
        pts = cond["xyz"].float().cpu()
        z = torch.stack([taubin_smooth(z[i], pts[i], torch.as_tensor(bnd[i][0], dtype=torch.float32),
                                       torch.as_tensor(bnd[i][1], dtype=torch.float32)) for i in range(B)])
    pe = sit.cond_embed(cond)

    def f(x, t):
        return sit(x, torch.full((B,), float(t), device=dev), cond, pe=pe)

    x, nfe = (dopri5(f, z.to(dev), atol=atol, rtol=rtol) if method == "dopri5"
              else fixed_step(f, z.to(dev), steps, method))
    return sit.denormalize(x), nfe


# =============================================================================
# 3. EVALUATION
# =============================================================================


def stack_cond(conds):
    return {k: torch.cat([c[k] for c in conds]) for k in conds[0]}


def cond_to(cond, dev):
    return {k: v.to(dev) for k, v in cond.items()}


def layout_sample(pack, i, vcfg, seed):
    """Layout i of the pack, unaugmented, as one batch of the SQ-VAE's encoder input."""
    return {k: v[None] for k, v in make_sample(pack.record(i), dict(vcfg, n_query=1), np.random.default_rng(seed),
                                               False).items()}


class LossEval:
    """The flow loss of a layout set at fixed draws: per layout, its unaugmented posterior mean, and k pairs (t, x0)
    with t stratified over (0, 1) and x0 from a seeded generator. Deterministic, so evaluations at different steps
    compare. Also reported by t bin, since the loss near t = 1 (little noise left) is where a sample's details are
    decided."""

    def __init__(self, pack, indices, vae, geom, vcfg, dev, k=4):
        self.items, self.k = [], k
        for n, i in enumerate(indices):
            x1, cond = encode_batch(vae, geom, tr.to_dev(layout_sample(pack, i, vcfg, 1000 + n), dev))
            x0 = torch.randn((k,) + x1.shape[1:], generator=torch.Generator().manual_seed(3000 + n))
            self.items.append((x1.cpu(), {k_: v.cpu() for k_, v in cond.items()}, x0))
        self.t = (torch.arange(k) + 0.5) / k

    @torch.no_grad()
    def __call__(self, sit, dev, chunk=2):
        was = sit.training
        sit.eval()
        per = []
        for s in range(0, len(self.items), chunk):
            part = self.items[s:s + chunk]
            x1 = sit.normalize(torch.cat([x.repeat(self.k, 1, 1) for x, _, _ in part]).to(dev))
            cond = cond_to(stack_cond([{k_: v.repeat(self.k, 1, 1) for k_, v in c.items()} for _, c, _ in part]), dev)
            x0 = torch.cat([x0 for _, _, x0 in part]).to(dev)
            t = self.t.repeat(len(part)).to(dev)
            per.append(flow_loss(sit, x1, cond, t, x0).cpu().reshape(len(part), self.k))
        sit.train(was)
        per = torch.cat(per)
        out = {"loss": per.mean().item()}
        out.update({f"loss_t{float(t):.2f}": per[:, j].mean().item() for j, t in enumerate(self.t)})
        return out


def wall_iou(pred, true):
    return float((pred & true).sum() / max(1, (pred | true).sum()))


def fps_cache(gcfg):
    """The Geom-AE run's token-set cache (domains.load_domains), or None."""
    path = os.path.join(from_repo(gcfg["out"]), f"fps_cache_N{gcfg['num_latents']}.npz")
    return np.load(path) if os.path.exists(path) else None


def domain_condition(geom, gcfg, obj, dev, seed=0, cache=None):
    """Sec. 4 steps 1-3 for one domain OBJ, as geom_ae.encode_obj does them: normalise, FPS N and 4N area samples,
    sdf there, then the Geom-AE code from its posterior mean. Returns the SiT condition on `dev`, the Domain (its
    frame: x_obj = x / scale + centre), and its boundary segments in that frame for the Taubin noise. The token sets
    are the first of the ones the Geom-AE trained with when `cache` (fps_cache()) holds the domain, which every corpus
    domain is; otherwise one set of each size is drawn. encode_obj draws five and keeps one, which cost 2 s a domain
    and made the evaluation sets take minutes to build."""
    d = Domain(obj)
    n = gcfg["num_latents"]
    key = os.path.basename(obj)
    if cache is not None and f"{key}/nf" in cache.files and int(cache[f"{key}/nf"]) == len(d.F):
        d.fps = {m: (cache[f"{key}/tri{m}"][:1], cache[f"{key}/bary{m}"][:1]) for m in (n, 4 * n)}
    else:
        rng = np.random.default_rng(seed)
        d.fps = {}
        for m in (n, 4 * n):
            tri, bary = geo.sample_area(d.V, d.F, d.cdf, max(8 * m, 8192), rng)
            sel = geo.farthest_point_sampling(geo.point_from_bary(d.V, d.F, tri, bary), m, rng)
            d.fps[m] = (tri[sel][None], bary[sel][None].astype(np.float32))
    s = geom_sample(d, dict(gcfg, n_query=1, n_near=1, n_global=1), np.random.default_rng(seed), False)
    b = tr.to_dev({k: v[None] for k, v in s.items() if "fps" in k}, dev)
    cond = {"code": geom_code(geom, b), "xyz": b[f"xyz_fps_{n}"], "sdf": b[f"sdf_fps_{n}"]}
    return cond, d, (d.V[d.E[:, 0]], d.V[d.E[:, 1]])


def to_frame(P, src, dst):
    """Points from one normalised frame (anything with centre / scale) to another."""
    return (P / src.scale + src.centre - dst.centre) * dst.scale


class DomainEval:
    """Samples per domain, scored against the closest of the domain's layouts (wall IoU at EVAL_POINTS area samples of
    the OBJ), with the SQ-VAE's reconstruction of each layout, scored the same way, as the ceiling. The condition
    comes from the OBJ as at inference."""

    def __init__(self, pack, indices, vae, geom, vcfg, gcfg, cfg, dev):
        by = {}
        for i in indices:
            by.setdefault(pack.rows[i]["model"], []).append(i)
        self.wall, self.domains = cfg["wall"], []
        cache = fps_cache(gcfg)
        for n, (name, idx) in enumerate(by.items()):
            rs = [pack.rows[i] for i in idx]
            cond, d, bnd = domain_condition(geom, gcfg, rs[0]["path"], dev, seed=n, cache=cache)
            rng = np.random.default_rng(2000 + n)
            P = geo.point_from_bary(d.V, d.F, *geo.sample_area(d.V, d.F, d.cdf, EVAL_POINTS, rng))
            walls, recon = [], []
            for i, r in zip(idx, rs):
                L = Layout(r["npz"])
                PL = to_frame(P, d, L)
                tri, bary, _ = locate(L.V, L.T, PL)
                walls.append(layout_fields(L, L.V, tri, bary)["cdf"][:, 0] <= self.wall)
                b = tr.to_dev(layout_sample(pack, i, vcfg, 0), dev)
                z = vae.encode(b, geom_code(geom, b), is_mode=True)["latent"]
                recon.append(decode_fields(vae, z, PL)["cdf"][:, 0] <= self.wall)
            ceiling = float(np.mean([max(wall_iou(p, w) for w in walls) for p in recon]))
            self.domains.append({"name": name, "rows": rs, "cond": {k: v.cpu() for k, v in cond.items()}, "domain": d,
                                 "bnd": bnd, "P": P.astype(np.float32), "walls": walls, "ceiling": ceiling})

    def score(self, vae, z, D, dev):
        """(wall IoU with the closest layout, its index) of latent z (N, C) on domain D."""
        pred = decode_fields(vae, z[None].to(dev), D["P"])["cdf"][:, 0] <= self.wall
        ious = [wall_iou(pred, w) for w in D["walls"]]
        return max(ious), int(np.argmax(ious))

    @torch.no_grad()
    def __call__(self, sit, vae, dev, n_gen=1, noise="plain", method="heun", steps=12, atol=1e-6, rtol=1e-3,
                 chunk=8, keep=False):
        was = sit.training
        sit.eval()
        per, nfes = [], []
        for s in range(0, len(self.domains), chunk):
            part = self.domains[s:s + chunk]
            cond = cond_to(stack_cond([D["cond"] for D in part]), dev)
            res = [[] for _ in part]
            for g in range(n_gen):
                seeds = [1_000_003 * g + s + j for j in range(len(part))]
                z, nfe = sample_latents(sit, cond, seeds, noise, [D["bnd"] for D in part], method, steps, atol, rtol)
                nfes.append(nfe)
                for j, D in enumerate(part):
                    iou, which = self.score(vae, z[j], D, dev)
                    res[j].append({"iou": iou, "layout": which, "z": z[j].cpu() if keep else None})
            for D, rs in zip(part, res):
                per.append({"name": D["name"], "samples": rs, "mean": float(np.mean([r["iou"] for r in rs])),
                            "best": max(r["iou"] for r in rs), "ceiling": D["ceiling"]})
        sit.train(was)
        out = {"iou": float(np.mean([p["mean"] for p in per])), "best": float(np.mean([p["best"] for p in per])),
               "ceiling": float(np.mean([p["ceiling"] for p in per])), "nfe": float(np.mean(nfes))}
        return (out, per) if keep else out


def plot_samples(ev, per, vae, dev, out, title, k=6, h_rel=0.004):
    """A row per domain, the worst best-of-n first, so the failures are the first thing seen. The first column is the
    true layout the best sample came closest to, with its chart sides; the others are the samples. Every panel is the
    CDF on the domain's OBJ refined to about h_rel x its diagonal, coloured per triangle at its centroid, in OBJ
    coordinates."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap

    order = list(np.argsort([p["best"] for p in per]))
    idx = order[:k // 2] + order[::-1][:k - k // 2][::-1] if len(order) > k else order
    n_gen = len(per[0]["samples"])
    meshes = [extraction_mesh(ev.domains[i]["rows"][0]["path"], h_rel)[:2] for i in idx]
    aspect = [float(np.clip(np.ptp(V[:, 1]) / np.ptp(V[:, 0]), 0.45, 1.3)) for V, _ in meshes]
    cmap = LinearSegmentedColormap.from_list("field", FIELD_RAMP)
    fig, axes = plt.subplots(len(idx), 1 + n_gen, figsize=(2.6 * (1 + n_gen), sum(2.3 * a + 0.65 for a in aspect) + 0.5),
                             squeeze=False, gridspec_kw={"height_ratios": [a + 0.28 for a in aspect]})
    im = None
    for row, i, (V, T) in zip(axes, idx, meshes):
        D, p = ev.domains[i], per[i]
        C = V[T].mean(1)
        best = max(p["samples"], key=lambda s: s["iou"])
        r = D["rows"][best["layout"]]
        L = Layout(r["npz"])
        tri, bary, _ = locate(L.V, L.T, (C - L.centre) * L.scale)
        fields = [layout_fields(L, L.V, tri, bary)["cdf"][:, 0]]
        Pd = ((C - D["domain"].centre) * D["domain"].scale).astype(np.float32)
        fields += [decode_fields(vae, s["z"][None].to(dev), Pd)["cdf"][:, 0] for s in p["samples"]]
        heads = [f"{D['name']} ({len(D['rows'])} layout" + ("" if len(D["rows"]) == 1 else "s")
                 + f", SQ-VAE IoU {p['ceiling']:.2f})\nclosest: {r['method']}, {r['charts']} chart"
                 + ("" if r["charts"] == "1" else "s")]
        heads += [f"sample {j + 1}\nIoU {s['iou']:.2f} vs {D['rows'][s['layout']]['method']}"
                  for j, s in enumerate(p["samples"])]
        for j, (ax, vals) in enumerate(zip(row, fields)):
            im = ax.tripcolor(V[:, 0], V[:, 1], T, facecolors=np.clip(vals, 0, 1), cmap=cmap, vmin=0, vmax=1,
                              shading="flat", edgecolors="none", rasterized=True)
            if j == 0:
                ax.add_collection(LineCollection(chart_sides(r["npz"], L.V) / L.scale + L.centre,
                                                 colors="#1a1a19", linewidths=0.6))
            ax.set_title(heads[j], fontsize=7.5, color="#52514e")
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color("#c9c8c4")
    fig.suptitle(f"{title}\nIoU: decoded vs true walls (cdf <= {ev.wall}) of the closest layout; "
                 "lines: its chart sides", y=0.997, fontsize=9, color="#1a1a19")
    fig.tight_layout(rect=(0, 0, 0.9, 0.97))
    cax = fig.add_axes((0.915, 0.25, 0.012, 0.5))
    cb = fig.colorbar(im, cax=cax, ticks=[0, 0.25, 0.5, 0.75, 1])
    cb.set_label("CDF", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.savefig(out, dpi=110, facecolor="#fcfcfb")
    plt.close(fig)


# =============================================================================
# 4. TRAINING
# =============================================================================


def load_vae(path, dev):
    """SQ-VAE, its Geom-AE, both configs and the step, from sq_vae.py's checkpoint."""
    vae, geom, vcfg = load_sq_vae(path, dev)
    ck = torch.load(path, map_location="cpu", weights_only=False)
    return vae, geom, vcfg, ck["geom_config"], ck["step"]


def fingerprint(*modules):
    """A hash of the weights of `modules`: what identifies the frozen SQ-VAE and Geom-AE an SQ-Diffuse checkpoint was
    trained on, wherever their files have been copied to. The bytes are the fp32 CPU copies, so it is the same on
    every device."""
    h = hashlib.sha1()
    for m in modules:
        for k, v in sorted(m.state_dict().items()):
            h.update(k.encode())
            h.update(v.detach().to("cpu", torch.float32).contiguous().numpy().tobytes())
    return h.hexdigest()[:16]


def save(path, sit, cfg, step, metrics):
    torch.save({"model": {k: v.detach().cpu() for k, v in sit.state_dict().items()}, "config": portable(cfg),
                "step": step, "metrics": metrics}, path)


def load_sq_diffuse(path, device="cpu", vae_path=None):
    """SQ-Diffuse with the frozen SQ-VAE and Geom-AE it was trained on: returns (sit, vae, geom, cfg, vae_cfg,
    geom_cfg). The SQ-VAE is the one the checkpoint records (its path relative to the repository), or vae_path;
    either way its weights must hash to the recorded vae_fingerprint, since a different SQ-VAE decodes the same latent
    into a different layout."""
    ck = torch.load(path, map_location=device, weights_only=False)
    cfg = ck["config"]
    vae_path = vae_path or from_repo(cfg["vae"])
    if not os.path.exists(vae_path):
        raise SystemExit(f"{path} was trained on the SQ-VAE {cfg['vae']}, which is not at {vae_path}: copy that run "
                         "directory there (keep the layout under py/squadgen/runs), or pass its best.pt (--vae)")
    vae, geom, vcfg, gcfg, vstep = load_vae(vae_path, device)
    want = cfg.get("vae_fingerprint")
    if want is not None and fingerprint(vae, geom) != want:
        raise SystemExit(f"{vae_path} is not the SQ-VAE {path} was trained on (step {cfg['vae_step']}, weights "
                         f"{want}): its weights differ. Pass the right one with --vae.")
    if want is None and vstep != cfg["vae_step"]:
        print(f"warning: {path} was trained on SQ-VAE step {cfg['vae_step']}, loaded step {vstep}", flush=True)
    sit = build_model(cfg).to(device)
    sit.load_state_dict(ck["model"])
    return sit.eval(), vae, geom, cfg, vcfg, gcfg


def lr_at(step, cfg):
    return cfg["lr"] if cfg["const_lr"] else tr.lr_at(step, cfg)


def train(cfg, sit, vae, geom, vcfg, pack, train_idx, evals, log):
    dev = cfg["device"]
    amp = tr.amp_dtype(dev, cfg["amp"])
    opt = torch.optim.AdamW(sit.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    scaler = tr.grad_scaler(dev, amp)
    stream = dict(vcfg, n_query=1, seed=cfg["seed"], augment=cfg["augment"])
    loader = tr.make_loader(LayoutStream(pack, train_idx, stream), cfg)
    log(f"amp={amp}  workers={cfg['workers']}")
    best, t0, run = float("inf"), time.time(), {}
    for step, batch in enumerate(loader, 1):
        for g in opt.param_groups:
            g["lr"] = lr_at(step - 1, cfg)
        x1, cond = encode_batch(vae, geom, tr.to_dev(batch, dev))     # fp32: the targets must be exact
        with tr.autocast(dev, amp):
            loss = flow_loss(sit, sit.normalize(x1), cond).mean()
        gnorm = tr.train_step(sit, opt, scaler, loss, cfg["clip_grad"])
        tr.accumulate(run, {"loss": loss, "gnorm": gnorm})
        if step % cfg["log_every"] == 0 or step == 1:
            el = time.time() - t0
            log(f"step {step:6d}  lr {opt.param_groups[0]['lr']:.2e}  {tr.averages(run)}"
                + f"  ({el / step * 1000:.0f} ms/step, {el / 60:.1f} min{tr.memory_note(dev)})")
            run = {}
        if step % cfg["eval_every"] == 0 or step == cfg["steps"]:
            metrics, parts = {}, []
            for name, le, de in evals:
                m = le(sit, dev)
                m.update(de(sit, vae, dev, steps=cfg["eval_steps"]))
                metrics[name] = m
                parts.append(f"{name} loss {m['loss']:.4f} (t=0.88: {m['loss_t0.88']:.4f})  sample IoU {m['iou']:.3f}"
                             f" (SQ-VAE {m['ceiling']:.3f})")
            log(f"  eval {step:6d}  " + "  |  ".join(parts))
            save(os.path.join(cfg["out"], "last.pt"), sit, cfg, step, metrics)
            if metrics[evals[0][0]]["loss"] < best:
                best = metrics[evals[0][0]]["loss"]
                save(os.path.join(cfg["out"], "best.pt"), sit, cfg, step, metrics)
            # MPS keeps the evaluation's blocks cached after it, and training then needs its own on top: on 8 GB that
            # tipped the 2026-10-29 run into swap (stalled after its first evaluation). The evaluations also run in
            # chunks of the training batch's size (LossEval 2 layouts x 4 draws, DomainEval 8 domains) so their blocks
            # fit the same sizes.
            tr.empty_cache(dev)
        if step >= cfg["steps"]:
            break


def final_eval(sit, vae, evals, cfg, step, log):
    """best.pt with the release's sampler: n_gen samples per domain by dopri5, from Taubin-smoothed noise
    (infer_sqdiffuse.py) and from white noise (engine_sqdiffuse.evaluate)."""
    dev = cfg["device"]
    rows = []
    for name, _, de in evals:
        for noise in ("taubin", "plain"):
            t0 = time.time()
            m, per = de(sit, vae, dev, n_gen=cfg["n_gen"], noise=noise, method="dopri5", atol=cfg["atol"],
                        rtol=cfg["rtol"], keep=True)
            closest = {}
            for p, D in zip(per, de.domains):
                for s in p["samples"]:
                    meth = D["rows"][s["layout"]]["method"]
                    closest[meth] = closest.get(meth, 0) + 1
                    rows.append({"split": name, "model": p["name"], "noise": noise, "iou": round(s["iou"], 4),
                                 "closest_method": meth, "closest_charts": D["rows"][s["layout"]]["charts"],
                                 "layouts": len(D["rows"]), "sq_vae_ceiling": round(p["ceiling"], 4)})
            log(f"best.pt (step {step}) on {name}, {len(per)} domains, {noise} noise, dopri5 ({m['nfe']:.0f} NFE per "
                f"solve, {time.time() - t0:.0f} s): sample IoU {m['iou']:.3f}, best of {cfg['n_gen']} {m['best']:.3f}, "
                f"SQ-VAE ceiling {m['ceiling']:.3f}  (min {min(p['best'] for p in per):.3f}, median best "
                f"{np.median([p['best'] for p in per]):.3f}); closest layouts: "
                + ", ".join(f"{k} {v}" for k, v in sorted(closest.items(), key=lambda kv: -kv[1])))
            if noise == "taubin":
                plot_samples(de, per, vae, dev, os.path.join(cfg["out"], f"samples_{name}.png"),
                             f"SQ-Diffuse 2D, best.pt (step {step}), {name} domains: {cfg['n_gen']} samples each "
                             f"(dopri5, Taubin-smoothed noise), worst and best")
    with open(os.path.join(cfg["out"], "samples.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="squadgen2d sq-diffuse", description=__doc__.split("\n")[1].strip())
    ap.add_argument("--vae", help="sq_vae.py checkpoint (SQ-VAE + Geom-AE), frozen (default: the one --out's best.pt "
                                  f"was trained on with --eval_only, else {to_repo(DEFAULT_VAE)})")
    ap.add_argument("--latent_stats", help="target_mu / target_std (default: latent_stats.json beside --vae)")
    ap.add_argument("--data", default=DEFAULT_CDF, help="a cdf-data output directory (manifest.csv + npz/)")
    ap.add_argument("--pack", help="its field pack (default: DATA/pack_N<the Geom-AE's N>; see cdf-pack)")
    ap.add_argument("--out", default=DEFAULT_SQ_DIFFUSE)
    ap.add_argument("--methods", nargs="+", default=list(METHODS), choices=METHODS)
    ap.add_argument("--overfit", nargs="+", metavar="MODEL",
                    help="train and evaluate on these domains' layouts only (Sec. 6 step 6)")
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--batch", type=int, default=8)                  # release: 1 per GPU x 16 GPUs; 8 GB: 8
    ap.add_argument("--depth", type=int, default=12)                 # Sec. 3.5: 12 x 6 x 64
    ap.add_argument("--n_heads", type=int, default=6)
    ap.add_argument("--d_head", type=int, default=64)
    ap.add_argument("--mlp_ratio", type=float, default=4.0)
    ap.add_argument("--lr", type=float, default=1e-4)                # release
    ap.add_argument("--const_lr", action="store_true", help="the release's constant rate, no warm-up or decay")
    ap.add_argument("--min_lr", type=float, default=1e-6)
    ap.add_argument("--warmup_frac", type=float, default=0.03)
    ap.add_argument("--wd", type=float, default=0.01)                # torch AdamW default, as the release
    ap.add_argument("--clip_grad", type=float, default=1.0)
    ap.add_argument("--no_aug", action="store_true")
    ap.add_argument("--wall", type=float, default=0.1, help="CDF wall threshold of the sample IoU (QuadExtraction -t)")
    ap.add_argument("--eval_steps", type=int, default=12, help="Heun steps of the training-time sample IoU")
    ap.add_argument("--n_gen", type=int, default=3, help="samples per domain in the final evaluation")
    ap.add_argument("--atol", type=float, default=1e-6)              # release sample_ode
    ap.add_argument("--rtol", type=float, default=1e-3)
    tr.add_device_args(ap)
    ap.add_argument("--eval_only", action="store_true", help="evaluate and plot best.pt")
    args = ap.parse_args(argv)

    best = os.path.join(args.out, "best.pt")
    if args.vae is None:
        # Evaluating a run (perhaps trained elsewhere and copied here): the SQ-VAE it records.
        recorded = (torch.load(best, map_location="cpu", weights_only=False)["config"]["vae"]
                    if args.eval_only and os.path.exists(best) else None)
        args.vae = from_repo(recorded) if recorded else DEFAULT_VAE
    if not os.path.exists(args.vae):
        raise SystemExit(f"no SQ-VAE checkpoint {args.vae}: train one with  python -m squadgen2d sq-vae, or pass --vae")
    stats_path = args.latent_stats or os.path.join(os.path.dirname(args.vae), "latent_stats.json")
    if not os.path.exists(stats_path):
        raise SystemExit(f"no {stats_path}: write it with  python -m squadgen2d sq-vae --eval_only --out "
                         f"{os.path.dirname(args.vae)}")
    dev = tr.pick_device(args.device)
    vae, geom, vcfg, gcfg, vstep = load_vae(args.vae, dev)
    with open(stats_path) as f:
        stats = json.load(f)
    if stats.get("step") != vstep:
        raise SystemExit(f"{stats_path} is for SQ-VAE step {stats.get('step')}, {args.vae} is step {vstep}: "
                         "rerun sq-vae --eval_only")
    pack = open_pack(args.data, args.pack, vcfg["n_tokens"])
    cfg = dict(vars(args), pack=pack.path, device=dev, augment=not args.no_aug, vae=os.path.abspath(args.vae),
               vae_step=vstep, vae_fingerprint=fingerprint(vae, geom), in_channels=vcfg["latent_dim"],
               in_channels_cond=gcfg["dim"], qk_norm=1, target_mu=stats["target_mu"], target_std=stats["target_std"],
               workers=tr.default_workers(dev) if args.workers is None else args.workers)

    log = tr.RunLog(args.out)
    log.header(cfg, skip=("target_mu", "target_std"))
    for r in pack.rows:
        r["path"] = os.path.join(REPO, r["path"])
    keep = [i for i, r in enumerate(pack.rows) if r["method"] in args.methods]
    if args.overfit:
        idx = [i for i in keep if pack.rows[i]["model"] in args.overfit]
        if not idx:
            raise SystemExit(f"no layouts of {' '.join(args.overfit)} in {pack.path}")
        train_idx, sets = idx, [("overfit", idx, idx)]
        log(f"overfit: {len(idx)} layouts of {' '.join(sorted({pack.rows[i]['model'] for i in idx}))}")
    else:
        # The SQ-VAE's split, so that SQ-Diffuse validates on the same held-out families.
        if vcfg.get("val_families") is not None:
            val_fams, how = set(vcfg["val_families"]), "recorded in the SQ-VAE checkpoint"
        else:
            val_fams, how = geom_ae_val_families(gcfg)
        tr_i, va_i = split_layouts([pack.rows[i] for i in keep], val_fams)
        train_idx, val_idx = [keep[i] for i in tr_i], [keep[i] for i in va_i]
        log(f"split: the manifest's split column (the Geom-AE's held-out families: {how})" if val_fams is None
            else f"split: the held-out families, {how}")
        # The 48 training layouts sq_vae.py's train eval used (same seed, same draw), and for the sample IoU every
        # training layout of their domains.
        rng = np.random.default_rng(args.seed)
        drawn = [train_idx[i] for i in rng.choice(len(train_idx), min(48, len(train_idx)), replace=False)]
        names = {pack.rows[i]["model"] for i in drawn}
        sets = [("val", val_idx, val_idx), ("train", drawn, [i for i in train_idx if pack.rows[i]["model"] in names])]
        log(f"{len(keep)} layouts of {len({pack.rows[i]['model'] for i in keep})} domains: {len(train_idx)} train, "
            f"{len(val_idx)} val ({len({pack.rows[i]['model'] for i in val_idx})} domains)")

    t0 = time.time()
    evals = [(name, LossEval(pack, lidx, vae, geom, vcfg, dev), DomainEval(pack, didx, vae, geom, vcfg, gcfg, cfg, dev))
             for name, lidx, didx in sets]
    log(f"evaluation sets ready in {time.time() - t0:.0f} s: "
        + ", ".join(f"{n} {len(le.items)} layouts / {len(de.domains)} domains" for n, le, de in evals))

    torch.manual_seed(args.seed)
    sit = build_model(cfg).to(dev)
    log(f"device={dev}  SiT parameters={sum(p.numel() for p in sit.parameters()):,}  hidden {args.n_heads * args.d_head}"
        f"  (frozen: SQ-VAE {sum(p.numel() for p in vae.parameters()):,}, Geom-AE {sum(p.numel() for p in geom.parameters()):,})")
    if not args.eval_only:
        train(cfg, sit, vae, geom, vcfg, pack, train_idx, evals, log)
    sit, vae, geom, _, _, _ = load_sq_diffuse(best, dev, args.vae)
    step = torch.load(best, map_location="cpu", weights_only=False)["step"]
    final_eval(sit, vae, evals, cfg, step, log)
    log(f"plots in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
