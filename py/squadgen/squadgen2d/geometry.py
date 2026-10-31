"""
Planar geometry in numpy: a domain is a triangle mesh of Omega at z = 0.

OBJ reading, area sampling and farthest-point sampling, the truncated SDF (docs/squadgen_2d.md Sec. 2.7 item 3), the
normalised frame of Sec. 4 step 1, and the training augmentation (Sec. 2.7 item 1). Nothing here knows about
networks or layouts; domains.py builds the Geom-AE's domains from it and layout.py the SQ-VAE's layouts.
"""

import numpy as np


def read_obj(path):
    """Vertices (V, 2) and triangles (F, 3) of a CrossGen OBJ. Material tags are ignored: the Geom-AE sees only the
    domain."""
    V, Fc = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                V.append(line.split()[1:3])
            elif line.startswith("f "):
                Fc.append([int(t.split("/")[0]) - 1 for t in line.split()[1:4]])
    return np.asarray(V, np.float64), np.asarray(Fc, np.int64)


def cross2(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def perp(a):
    return np.stack([a[..., 1], -a[..., 0]], -1)


def signed_areas(V, Fc):
    a, b, c = V[Fc[:, 0]], V[Fc[:, 1]], V[Fc[:, 2]]
    return 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))


def boundary_edges(Fc):
    """Edges used by exactly one triangle, each as its triangle runs it. They close up into the loops of dOmega, so
    the even-odd test below needs no orientation, and on a consistently oriented mesh they run with Omega on one
    side, which boundary_frame() relies on."""
    E = np.concatenate([Fc[:, [0, 1]], Fc[:, [1, 2]], Fc[:, [2, 0]]])
    _, inv, cnt = np.unique(np.sort(E, 1), axis=0, return_inverse=True, return_counts=True)
    return E[cnt[inv.ravel()] == 1]


def normalisation(V, Fc, bverts):
    """The inference frame (Sec. 4 step 1) as (centre, scale): area centroid at the origin, farthest boundary vertex
    at radius 0.999. The area centroid rather than the vertex mean, so that the frame does not depend on how the
    domain was meshed."""
    ar = np.abs(signed_areas(V, Fc))
    c = (ar[:, None] * V[Fc].mean(1)).sum(0) / ar.sum()
    return c, 0.999 / np.linalg.norm(V[bverts] - c, axis=1).max()


def normalise(V, Fc, bverts):
    c, s = normalisation(V, Fc, bverts)
    return (V - c) * s


def boundary_frame(A, B):
    """normalisation() from dOmega alone, as oriented segments A -> B: the area centroid of the region they bound
    (the shoelace sums, holes included because they run the other way) and 0.999 over the largest distance of a
    segment start from it. Every vertex of dOmega starts one segment, so this is the frame normalisation() gives
    the triangulated domain, up to rounding."""
    w = cross2(A, B)
    c = ((A + B) * w[:, None]).sum(0) / (3 * w.sum())
    return c, 0.999 / np.linalg.norm(A - c, axis=1).max()


def point_from_bary(V, Fc, tri, bary):
    return np.einsum("pk,pkd->pd", bary, V[Fc[tri]])


def sample_area(V, Fc, cdf, n, rng):
    """n uniform samples of Omega as (triangle, barycentric) pairs. `cdf` is the cumulative triangle area; it is a
    ratio, so it survives the similarity part of any augmentation and the anisotropic part only changes the density
    by up to 10 %."""
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
    """Points per chunk so that a (chunk, E, 2) temporary stays near 16 MB. Several DataLoader workers each hold a
    few of them, on machines with 8 GB."""
    return max(16, (1 << 20) // max(1, n_edges))


def _point_segment(p, a, b):
    e = b - a
    t = np.clip(((p - a) * e).sum(-1) / np.maximum((e * e).sum(-1), 1e-30), 0.0, 1.0)
    return np.sqrt(((p - a - t[..., None] * e) ** 2).sum(-1))


def segment_distance(x, A, B, cap=np.inf):
    """Distance from each x (P, 2) to the segment set {[A_e, B_e]}: exact where it is below `cap`, and `cap`
    elsewhere. The truncated SDF needs nothing beyond tau, and a finite cap lets each point look only at the edges
    registered in its cell of a cap-sized grid, about a tenth of dOmega. That is what keeps two DataLoader workers
    ahead of the GPU. An edge is registered in every cell that its bounding box, grown by cap, touches. A point
    within cap of the edge lies in that grown box, so its cell has the edge."""
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
    """clip(phi / tau, -1, 1) as a (P, 1) column, phi < 0 inside Omega. Area samples are inside by construction, and
    skipping the ray cast for them keeps a sample that lands on dOmega from being misclassified by rounding."""
    d = segment_distance(x, A, B, cap=tau)
    d = -d if known_inside else np.where(inside(x, A, B), -d, d)
    return np.clip(d / tau, -1.0, 1.0).astype(np.float32)[:, None]


# =============================================================================
# AUGMENTATION  (Sec. 2.7 item 1 / Appendix D add_noise_2d, without the translation: the frame re-centres anyway)
# =============================================================================


class Deformation:
    """One random deformation of the plane: with probability 1/2 a smooth warp (three random low-frequency
    sinusoids, amplitude up to warp_amp; they play the part of the release's Perlin noise without the `noise`
    package), then a rotation about z, a reflection with probability 1/2, and an anisotropic scale in [0.9, 1.1]^2.

    It is a map of points, so it applies the same way to a mesh's vertices (augment()) and to points that carry
    fields with no mesh at all (cdf_pack.make_sample(), which also needs jacobian() to turn gradients). The draws
    are made in the order the mesh augmentation always made them, so one generator state gives the same
    deformation either way."""

    def __init__(self, rng, warp_amp=0.06):
        self.waves = []
        if rng.random() < 0.5:
            amp = rng.uniform(0, warp_amp)
            for _ in range(3):
                k = rng.normal(0, 1.5, 2)
                ph = rng.uniform(0, 2 * np.pi)
                self.waves.append((amp * rng.normal(0, 1, 2), k, ph))
        th = rng.uniform(-np.pi, np.pi)
        self.rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        self.reflect = rng.random() < 0.5
        self.scale = rng.uniform(0.9, 1.1, 2)

    def __call__(self, x):
        p = x.copy()
        for a, k, ph in self.waves:
            p += a * np.sin(x @ k + ph)[:, None] / 3
        p = p @ self.rot.T
        if self.reflect:
            p[:, 0] = -p[:, 0]
        return p * self.scale

    def jacobian(self, x):
        """d(self(x))/dx as (P, 2, 2): the warp's I + sum (a/3) k^T cos(k.x + ph), then the linear part."""
        J = np.broadcast_to(np.eye(2), (len(x), 2, 2)).copy()
        for a, k, ph in self.waves:
            J += np.cos(x @ k + ph)[:, None, None] * np.outer(a / 3, k)[None]
        L = self.rot.copy()
        if self.reflect:
            L[0] = -L[0]
        return (self.scale[:, None] * L)[None] @ J


def augment(V, Fc, rng, warp_amp=0.06):
    """A Deformation of a mesh's vertices. One that folds a triangle is rejected and another is drawn, up to five
    times, after which the mesh is left as it is. A reflection flips every triangle, which the test accepts."""
    a0 = signed_areas(V, Fc)
    for _ in range(5):
        p = Deformation(rng, warp_amp)(V)
        s = signed_areas(p, Fc) * a0
        if np.all(s > 0) or np.all(s < 0):
            return p
    return V.copy()
