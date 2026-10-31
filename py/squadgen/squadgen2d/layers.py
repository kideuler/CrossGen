"""
The network modules the Geom-AE and SQ-VAE share: the release's models_ae.py building blocks with the Sec. 3.1 edits
for 2D (POS_DIM = 2, every vector key 2-d). Module names and shapes are the release's, so that state-dict keys match
it; SQ-Diffuse's SiT reuses PointEmbed for its geometry conditioning.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

POS_DIM = 2             # Sec. 3.1 item 1
# Channels of every point feature and output key: the release's network/utils.py KEY_DIM_DICT with the Sec. 3.1
# item 1 edits. One table for both networks, as there: the Geom-AE reads the sdf keys, SQ-VAE the fields, their unit
# gradients and the seven offsets, every vector of them 2-d.
KEY_DIM_DICT = {"sdf": 1, "sdf_near": 1, "sdf_global": 1, "dcdf": 1, "cdf": 1, "gdcdf": 2, "gcdf": 2,
                "offset": 2, "offset1": 2, "offset2": 2, "offset3": 2, "offsetb": 2, "offsetc": 2, "offsetd": 2}


def xavier_init(module):
    """The release's initialize_weights base: Xavier-uniform weights and zero biases on every Linear."""
    def basic(m):
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)
    module.apply(basic)


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
    """Learned Fourier features of the position. mlp_pos holds random frequencies (initialised N(0, 1), times 2 pi)
    and is trained. One linear layer per point feature, and a linear layer over the concatenation. `query` is the
    position part alone; the decoder embeds its query points with it, so they share mlp_pos and mlp_pos_post with
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
