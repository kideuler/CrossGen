"""
A layout: one CDF NPZ's split quad mesh, and the chart distance fields on it (docs/squadgen_2d.md Appendix B).

cdf_data.py writes one NPZ per layout (a method's block decomposition of a domain, the blocks as charts; the release's
keys, written by BlockDecomposition.cdf() in C++). Layout reads its split mesh into the Geom-AE's frame, and
layout_fields() evaluates every field SQ-VAE reads at any points of it. cdf_pack.py evaluates them once at every
layout's token and query points for training; the evaluations and tools that need fields elsewhere (the round trip,
SQ-Diffuse's scoring, the plots) call layout_fields() directly.
"""

import csv
import os

import numpy as np

from . import geometry as geo

# The manifest's method names: CrossGen's five, and "oracle" (ORACLE's choice among them, Mesh.oracle()).
METHODS = ("zipline", "umber", "meridian", "torsion", "atlas", "oracle")
# The seven offsets in the order of the NPZ's corner ids: offset_abcd_id (own chart corner = dual-chart centre, side
# mid-point, chart centre, side mid-point), then offset_123_id (the chart's other three corners). Sec. 2.4.
OFFSET_KEYS = ("offset", "offsetb", "offsetc", "offsetd", "offset1", "offset2", "offset3")
NPZ_KEYS = ("subdiv_vertex", "subdiv_facet", "quad_split", "face2quad", "offset_abcd_id", "offset_123_id",
            "edge_vertex_color", "edge_color_store", "sample_bary_coords", "sample_point_tri_id")


def read_manifest(data, methods=METHODS):
    """The status-ok rows of DATA/manifest.csv for `methods`, in file order, with `npz` made absolute."""
    path = os.path.join(data, "manifest.csv")
    if not os.path.exists(path):
        raise SystemExit(f"no {path}: make it with  python -m squadgen2d cdf-data run --out {data}")
    with open(path) as f:
        rows = [r for r in csv.DictReader(f) if r["status"] == "ok" and r["method"] in methods]
    if not rows:
        raise SystemExit(f"no ok layouts of {' '.join(methods)} in {path}")
    for r in rows:
        r["npz"] = os.path.join(data, r["npz"])
    return rows


class Layout:
    """One NPZ's split mesh as triangles, with the vertices in the Geom-AE's frame. Each triangle carries what its
    fields need: the coordinate pairs of its two quad edges (Sec. 2.2 step 5) and the vertex ids of the seven offset
    targets. `centre`/`scale` map the frame back, x_obj = x / scale + centre, as domains.Domain does."""

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
        # A triangle's edges v0->v1 and v1->v2 are quad edges, the third is the diagonal. Each quad half-edge holds a
        # signed index into the (start, end) pairs; a negative index is the pair read backwards
        # (load_patch_data_reformat, Appendix B).
        e = np.tile([0, 1, 2, 3], (len(Q), 1))
        e[split] = [3, 0, 1, 2]
        e = e.reshape(-1, 2) + 4 * (np.arange(len(self.T)) // 2)[:, None]
        cid = a["edge_vertex_color"].ravel()[e]
        col = a["edge_color_store"].astype(np.float64)[np.abs(cid)]
        col[cid < 0] = col[cid < 0][:, ::-1]
        self.col = col
        sub = a["face2quad"].ravel()[np.arange(len(self.T)) // 2]
        self.corners = np.concatenate([a["offset_abcd_id"], a["offset_123_id"]], 1).astype(np.int64)[sub]
        self.E = geo.boundary_edges(self.T)
        self.bverts = np.unique(self.E)
        V = a["subdiv_vertex"][:, :2].astype(np.float64)
        self.centre, self.scale = geo.normalisation(V, self.T, self.bverts)
        self.V = (V - self.centre) * self.scale
        self.samples = (a["sample_point_tri_id"].ravel().astype(np.int64), a["sample_bary_coords"].astype(np.float64))


def layout_fields(L, V, tri, bary):
    """Every point feature SQ-VAE reads, at the points (tri, bary) of layout L with vertices V. dcdf, gdcdf, cdf and
    gcdf are Appendix B's: the release's get_patch_grading_colors with scalar cross products, equal to it to float32
    rounding. The two subchart coordinates c0, c1 are oblique projections onto the triangle's two quad edges, mapped
    to their coordinate intervals; then CDF = 1 - max(c0, c1) and DCDF = min(c0, c1). The seven offsets are corner
    minus point. Reflections need no care: both cross products in t0, t1 change sign together, and the gradients
    turn with the mesh."""
    t = L.T[tri]
    v0, v1, v2 = V[t[:, 0]], V[t[:, 1]], V[t[:, 2]]
    p = bary[:, 0:1] * v0 + bary[:, 1:2] * v1 + bary[:, 2:3] * v2
    e0, e1 = v1 - v0, v2 - v1
    den = geo.cross2(e0, e1)
    col = L.col[tri]
    t0 = np.clip(geo.cross2(p - v0, e1) / den, 0, 1)            # along e0, parallel to e1
    t1 = np.clip(geo.cross2(p - v1, e0) / -den, 0, 1)           # along e1, parallel to e0
    c0 = col[:, 0, 0] + t0 * (col[:, 0, 1] - col[:, 0, 0])
    c1 = col[:, 1, 0] + t1 * (col[:, 1, 1] - col[:, 1, 0])
    g0 = ((col[:, 0, 1] - col[:, 0, 0]) / den)[:, None] * geo.perp(e1)
    g1 = ((col[:, 1, 1] - col[:, 1, 0]) / -den)[:, None] * geo.perp(e0)
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


def chart_sides(path, V):
    """Segments (S, 2, 2) of the true chart sides: split-mesh edges on dOmega or between faces of two charts
    (subchart_block, which BlockDecomposition.cdf() writes)."""
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
