"""DeepCAD sketch profiles as gmsh .geo domains: more single-material models for data/meshes.

DeepCAD (Wu, Xiao & Zheng, ICCV 2021; the code in ~/Documents/DeepCAD, its data in
~/Documents/data) is Onshape models, one JSON per model under cad_json/<4 digits>/
<8 digits>.json, each a sequence of sketches and extrusions: 215,093 files, of which the
178,238 in train_val_test_split.json are the models of DeepCAD's paper (the others are in
no split, `none` in the manifest, and are converted all the same). DeepCAD's own reader
(cadlib/extrude.py: CADSequence.from_dict -> Extrude.from_dict -> sketch.Profile.from_dict)
takes, for every ExtrudeFeature, the sketch profiles it extrudes. A profile is one closed
region of its sketch's plane: loops of Line3D, Arc3D and Circle3D curves in the sketch's
own 2D frame, one loop around the outside and one per hole. cadlib/visualize.py
create_profile_face makes each profile an OpenCASCADE planar face and prisms it into the
model's solid, and dataset/json2pc.py samples point clouds off those solids. The 2D
geometry in DeepCAD is therefore the profile, and this script writes one .geo per profile:
a single-material domain like those in data/geometry/singlemat and MAMBO's planar faces.

    python3 py/deepcad_to_geo.py                      # every model -> data/geometry/deepcad
    python3 py/deepcad_to_geo.py --limit 2000 --out /tmp/dc   # the first 2000 models only

### What DeepCAD's reader gets wrong or leaves loose

* `is_outer` is true on every loop (all 17,937 loops of a 3000-model sample), so it does not
  say which loop is the outside. Here the outer loop is the one enclosing the most area, as
  in SKELETON's Step 1.
* An Arc3D runs counter-clockwise in the sketch's own frame from reference_vector, by
  end_angle (start_angle is always 0). That is how onshape-cad-parser (parser.py,
  _parse_edge_to_fusion360_format), which wrote these files, builds it: it sorts the two
  ends' angles, takes the counter-clockwise arc between them, and switches to the other arc
  when Onshape's midpoint of the edge is not at the average angle. DeepCAD's
  Arc.get_mid_point reads it the same way. The `normal` does not give the direction (a tenth
  of the arcs have a normal against the sketch's z axis and still run counter-clockwise in
  its frame), and start_point and end_point are the two ends in either order (a third are
  swapped). So the point half way round decides which of the two arcs between the ends is
  meant; the ends alone cannot, and a semicircle has two.
  The parser compares the midpoint angles rounded to three decimals, so it can switch an arc
  that should not have been switched (model 00005116 has a semicircle on the wrong side). A
  profile that this breaks is no longer a valid region, and the checks below reject it.
* A loop's curves are in chain order, but each may be reversed (sketch.Loop.reorder fixes
  this by comparing neighbours). Here ends within 1e-6 of the profile's size, or 1e-7 (in
  the files' metres; the parser rounds coordinates to 1e-8), are welded into one point, and
  each loop is walked as a cycle of those points.

By default the profiles written are those an ExtrudeFeature uses, the regions DeepCAD
builds its models from; --all-profiles adds the sketch regions no feature extruded. A
sketch's regions tile what it encloses, so most of those are copies of a written shape
anyway (the disk inside an extruded annulus).

### The .geo

It uses the OpenCASCADE kernel, as geom012.geo does. It holds Points, Lines, circle arcs
Circle(n) = {start, centre, end}, whole circles Circle(n) = {x, y, 0, r}, one Curve Loop per
loop and one Plane Surface, and it deletes the arc centres last. The centres have to go
because Mesh2Dgmsh normalises by gmsh's bounding box of every entity, points included, and
the centre of a shallow concave arc can lie far outside the face. That would shrink the face
inside [-1,1]^2 and coarsen its mesh. The built-in kernel cannot delete them, because its
arcs keep using their centres, so OCC is used. An arc of more than 135 degrees is written
as equal pieces of at most 135 degrees each. Both kernels read {start, centre, end} as the
shorter arc, which cannot express a semicircle or anything longer; the joins between the
pieces are smooth, so no corner is added.

Each profile is posed as Mesh3Dto2Dfacesgmsh poses a MAMBO face: its minimum-area bounding
rectangle along the axes, long side on x spanning [-1, 1]. Mesh2Dgmsh's own normalisation is
then the identity on it, and two copies of one shape get one mesh size, whichever way
round they were sketched.

### Which profiles are written

The rules are Mesh3Dto2Dfacesgmsh's, which picked the MAMBO faces. They are read off the
boundary as gmsh will mesh it at h = 2/--np: every curve is cut into gmsh's count of equal
pieces, max(1, floor(L/h + 0.99)), and at least round(7 sweep / 2 pi) on an arc
(Mesh.MinimumCircleNodes). Both were measured on gmsh 4 with Mesh2Dgmsh's uniform size. On
that polyline BoundaryFeatures' rules apply as written (a corner turns by more than 20
degrees; a run goes from one corner to the next, or is a whole loop with at most one corner).
A profile is written when it is

* a valid region: each loop is one closed cycle, the holes are inside the outer loop and
  outside one another, no two loops share a point, and the boundary polyline does not cross
  itself (if it did, gmsh's boundary mesh would cross too);
* nice: its shortest corner-to-corner run is at least --min-run mesh edges (3) and its
  thickness 2A/P at least --min-width mesh edges (4), as for MAMBO; and, beyond MAMBO's
  rules, no cusp, a corner sharper than 1 degree. There the boundary meets itself
  tangentially: an arc running into a line it touches (2.4 % of all extruded profiles), or
  the parser's wrongly switched fillet, which doubles back. No block layout has a block for
  a corner of no angle.

The whole set (2026-10-30, 13 min on 8 processes): 788,324 extruded profiles gave 97,784
.geo files (400 MB); 591,143 were copies of a shape already written (one disk shape alone
151,344 times), 98,254 not nice (55,275 a short feature, 24,144 thin, 18,835 a cusp) and
1,143 invalid (659 of them a hole outside the outer loop, mostly the parser's switched
arcs). 30,033 of the domains have at most four curves, 9,183 of them a single four-sided
loop.
* new: no profile written before it is the same shape, with Mesh3Dto2Dfacesgmsh's signature
  and sameShape tolerances (loops; perimeter and the two principal second moments of area,
  made scale-free by the area; sorted corner angles, run lengths and run turnings). Models
  go in id order, so the copy kept is the one in the lowest-numbered model, and `copies`
  in the manifest counts how many profiles were that shape.

### Output

OUT/<model id>_<k>.geo, where k is the profile's index among all of the model's profiles
(its sketches in sequence order, each sketch's profiles in the file's order), whether or not
the others are written. The name's family (cdf_data.family: everything before the last "_")
is then the DeepCAD model, so one model's faces fall on one side of the Geom-AE / SQ-VAE
family split. OUT/manifest.csv has a row for each .geo written. A rerun first deletes the
.geo files in OUT that this script wrote (their first line starts with "// DeepCAD ") and
leaves any other file alone.

make_meshes.sh meshes every data/geometry/<kind>/*.geo into data/meshes/<kind>, so the
default OUT ends up in data/meshes/deepcad. To mesh only these, in parallel:

    ls data/geometry/deepcad/*.geo | xargs -P 8 -I{} build/Mesh2Dgmsh {} data/meshes/deepcad 100

(Mesh2Dgmsh's "Unknown GEO point" error on these files is its built-in-kernel normalisation
failing over to the OCC one, as on geom012.geo; it is harmless.)
"""
import argparse
import bisect
import csv
import glob
import json
import math
import os
import sys
import time
from collections import Counter
from multiprocessing import Pool

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_DATA = os.path.expanduser("~/Documents/data")
DEFAULT_OUT = os.path.join(ROOT, "data", "geometry", "deepcad")
HEADER = "// DeepCAD "             # first line of every .geo this writes (see the docstring's Output)

CORNER = math.radians(20.0)        # BoundaryFeatures::Options::cornerAngle
MAX_ARC = math.radians(135.0)      # longest arc piece written, below the pi the kernels allow
WELD = 1e-6                        # ends this close (x the profile's bounding-box diagonal) are one point,
WELD_ABS = 1e-7                    # or this close in the files' metres (the parser rounds to 1e-8)
ARC_END = 1e-3                     # an arc's ends off its angles by more than this x radius: unreadable
MIN_CIRCLE_NODES = 7               # gmsh's Mesh.MinimumCircleNodes default
CUSP = math.radians(179.0)         # a turn this sharp leaves a corner of under 1 degree


class Reject(Exception):
    """A profile not written: kind is 'invalid' or 'not nice', why says what."""

    def __init__(self, kind, why):
        super().__init__(why)
        self.kind, self.why = kind, why


# ----------------------------------------------------------------------------------------
# Reading a profile
# ----------------------------------------------------------------------------------------

def xy(p):
    return np.array([p["x"], p["y"]], dtype=float)


def turned(v, a):
    """v turned counter-clockwise by a."""
    c, s = math.cos(a), math.sin(a)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])


def wrap(a):
    """a in [0, 2 pi)."""
    return a % (2.0 * math.pi)


def read_curves(loops):
    """The profile's loops as raw curves in its sketch's frame:
    ('L', start, end), ('A', start, end, centre, mid) or ('C', centre, r)."""
    out = []
    for loop in loops:
        curves = []
        for c in loop["profile_curves"]:
            kind = c.get("type")
            if kind == "Line3D":
                curves.append(("L", xy(c["start_point"]), xy(c["end_point"])))
            elif kind == "Circle3D":
                r = float(c["radius"])
                if not r > 0.0:
                    raise Reject("invalid", "a circle of radius %g" % r)
                curves.append(("C", xy(c["center_point"]), r))
            elif kind == "Arc3D":
                ctr, r, ref = xy(c["center_point"]), float(c["radius"]), xy(c["reference_vector"])
                a0, a1 = float(c["start_angle"]), float(c["end_angle"])
                s, e = xy(c["start_point"]), xy(c["end_point"])
                if not r > 0.0 or not np.linalg.norm(ref) > 0.0:
                    raise Reject("invalid", "a degenerate arc")
                ref = ref / np.linalg.norm(ref)
                p0, p1 = ctr + r * turned(ref, a0), ctr + r * turned(ref, a1)
                off = min(max(np.linalg.norm(p0 - s), np.linalg.norm(p1 - e)),
                          max(np.linalg.norm(p0 - e), np.linalg.norm(p1 - s)))
                if off > ARC_END * r:
                    raise Reject("invalid", "an arc whose ends are not at its angles")
                curves.append(("A", s, e, ctr, ctr + r * turned(ref, 0.5 * (a0 + a1))))   # its mid point
            else:
                raise Reject("invalid", "a %s curve" % kind)
        if curves:
            out.append(curves)
    if not out:
        raise Reject("invalid", "no curves")
    return out


class Region:
    """One profile, welded and walked: points (n x 2) and loops, each a list of directed
    edges ('L', i, j), ('A', i, j, centre, ccw) or ('C', centre, r, ccw), the region on the
    left (outer loop counter-clockwise, holes clockwise)."""

    def __init__(self, points, loops):
        self.points, self.loops = points, loops

    def transformed(self, f, scale):
        """The region under p -> f(p), a rotation-and-scale by `scale` (no mirror)."""
        def edge(e):
            if e[0] == "L":
                return e
            if e[0] == "A":
                return ("A", e[1], e[2], f(e[3]), e[4])
            return ("C", f(e[1]), e[2] * scale, e[3])
        return Region(np.array([f(p) for p in self.points]), [[edge(e) for e in L] for L in self.loops])


def reversed_loop(loop):
    out = []
    for e in reversed(loop):
        if e[0] == "L":
            out.append(("L", e[2], e[1]))
        elif e[0] == "A":
            out.append(("A", e[2], e[1], e[3], not e[4]))
        else:
            out.append(("C", e[1], e[2], not e[3]))
    return out


def weld(raw):
    """Points for the curves' ends, ends within WELD of the profile's size made one."""
    ends, extent = [], []
    for curves in raw:
        for c in curves:
            if c[0] == "C":
                extent += [c[1] - c[2], c[1] + c[2]]
            else:
                ends += [c[1], c[2]]
                extent += [c[1], c[2]]
    extent = np.array(extent)
    diag = float(np.linalg.norm(extent.max(0) - extent.min(0)))
    if not diag > 0.0:
        raise Reject("invalid", "no extent")
    if not ends:
        return np.zeros((0, 2)), [], diag
    E = np.array(ends)
    tol = max(WELD * diag, WELD_ABS)
    D = np.linalg.norm(E[:, None, :] - E[None, :, :], axis=2)
    label = -np.ones(len(E), dtype=int)
    points = []
    for i in range(len(E)):
        if label[i] >= 0:
            continue
        near = (D[i] <= tol) & (label < 0)
        label[near] = len(points)
        points.append(E[near].mean(0))
    return np.array(points), list(label), diag


def arc_edges(points, i, j, ctr, mid):
    """The arc i -> j through `mid`, with its centre moved onto the bisector of the welded
    ends (so that both are exactly one radius from it: OCC builds the circle through the
    start and finds the end on it), cut into pieces of at most MAX_ARC. Appends the cut
    points to `points`; returns the edges."""
    p, q = points[i], points[j]
    m, d = 0.5 * (p + q), q - p
    n = np.array([-d[1], d[0]]) / np.linalg.norm(d)
    c = m + n * float(np.dot(ctr - m, n))
    ang = lambda v: math.atan2(v[1] - c[1], v[0] - c[0])
    tp, tq, tm = ang(p), ang(q), ang(mid)
    ccw = wrap(tm - tp) < wrap(tq - tp)
    sweep = wrap(tq - tp) if ccw else wrap(tp - tq)
    pieces = max(1, math.ceil(sweep / MAX_ARC - 1e-9))
    r = float(np.linalg.norm(p - c))
    ids = [i]
    for k in range(1, pieces):
        t = tp + (1 if ccw else -1) * sweep * k / pieces
        points.append(c + r * np.array([math.cos(t), math.sin(t)]))
        ids.append(len(points) - 1)
    ids.append(j)
    return [("A", ids[k], ids[k + 1], c, ccw) for k in range(pieces)]


def build_region(raw):
    """Weld the ends, walk each loop as one cycle, pick the outer loop and orient them."""
    pts, label, diag = weld(raw)
    points = [p for p in pts]
    loops, at = [], 0
    for curves in raw:
        if len(curves) == 1 and curves[0][0] == "C":
            loops.append([("C", curves[0][1], curves[0][2], True)])
            continue
        if any(c[0] == "C" for c in curves):
            raise Reject("invalid", "a circle in a loop of several curves")
        edges = []
        for c in curves:
            i, j = label[at], label[at + 1]
            at += 2
            if i == j:
                # A zero-length line (DeepCAD's create_edge_3d drops them too), or an arc
                # welded shut, which a profile only has as a sliver.
                if c[0] == "A" and np.linalg.norm(c[4] - pts[i]) > 1e3 * max(WELD * diag, WELD_ABS):
                    raise Reject("invalid", "an arc that closes on itself")
                continue
            edges.append((c, i, j))
        if not edges:
            raise Reject("invalid", "a loop of zero length")
        # Each point of a cycle starts one edge and ends another.
        incident = {}
        for k, (_, i, j) in enumerate(edges):
            incident.setdefault(i, []).append(k)
            incident.setdefault(j, []).append(k)
        if any(len(v) != 2 for v in incident.values()):
            raise Reject("invalid", "a loop that is not one closed cycle")
        walk, used, v = [], set(), edges[0][1]
        k = 0
        while k not in used:
            used.add(k)
            c, i, j = edges[k]
            if i != v:
                i, j = j, i
                c = (c[0], c[2], c[1]) + c[3:]
            walk.append((c, i, j))
            v = j
            k = next(x for x in incident[v] if x != k) if len(edges) > 1 else k
        if len(used) != len(edges) or v != walk[0][1]:
            raise Reject("invalid", "a loop that is not one closed cycle")
        loop = []
        for c, i, j in walk:
            loop += [("L", i, j)] if c[0] == "L" else arc_edges(points, i, j, c[3], c[4])
        loops.append(loop)
    region = Region(np.array(points), loops)
    shared = Counter(v for L in loops for v in {e[k] for e in L if e[0] != "C" for k in (1, 2)})
    if any(n > 1 for n in shared.values()):
        raise Reject("invalid", "boundary loops touch")
    # The outer loop encloses the most area; then outer counter-clockwise, holes clockwise.
    areas = [signed_area(fine_samples(region, L)) for L in loops]
    outer = int(np.argmax(np.abs(areas)))
    oriented = []
    for k in [outer] + [k for k in range(len(loops)) if k != outer]:
        want = 1.0 if k == outer else -1.0
        oriented.append(loops[k] if areas[k] * want > 0 else reversed_loop(loops[k]))
    region.loops = oriented
    return region


# ----------------------------------------------------------------------------------------
# Sampling, pose
# ----------------------------------------------------------------------------------------

def arc_angles(region, e):
    c = e[3]
    p, q = region.points[e[1]], region.points[e[2]]
    t0 = math.atan2(p[1] - c[1], p[0] - c[0])
    t1 = math.atan2(q[1] - c[1], q[0] - c[0])
    sweep = wrap(t1 - t0) if e[4] else wrap(t0 - t1)
    return c, float(np.linalg.norm(p - c)), t0, (1.0 if e[4] else -1.0) * sweep


def edge_samples(region, e, pieces):
    """`pieces` equal steps along edge e, from its start (included) to its end (not)."""
    t = np.arange(pieces) / pieces
    if e[0] == "L":
        p, q = region.points[e[1]], region.points[e[2]]
        return p[None, :] + t[:, None] * (q - p)[None, :]
    if e[0] == "A":
        c, r, t0, dt = arc_angles(region, e)
    else:
        c, r, t0, dt = e[1], e[2], 0.0, (2.0 if e[3] else -2.0) * math.pi
    a = t0 + dt * t
    return c[None, :] + r * np.stack([np.cos(a), np.sin(a)], axis=1)


def edge_length(region, e):
    if e[0] == "L":
        return float(np.linalg.norm(region.points[e[2]] - region.points[e[1]]))
    if e[0] == "A":
        _, r, _, dt = arc_angles(region, e)
        return r * abs(dt)
    return 2.0 * math.pi * e[2]


def fine_samples(region, loop):
    """A loop finely enough for areas and the pose: lines as they are, arcs every degree."""
    out = []
    for e in loop:
        if e[0] == "L":
            out.append(edge_samples(region, e, 1))
        else:
            sweep = 2.0 * math.pi if e[0] == "C" else abs(arc_angles(region, e)[3])
            out.append(edge_samples(region, e, max(2, math.ceil(math.degrees(sweep)))))
    return np.concatenate(out)


def mesh_samples(region, loop, h):
    """A loop as gmsh meshes it at size h: each curve in max(1, floor(L/h + 0.99)) equal
    pieces, and an arc in at least round(7 sweep / 2 pi) (measured; see the docstring)."""
    out = []
    for e in loop:
        n = max(1, math.floor(edge_length(region, e) / h + 0.99))
        if e[0] != "L":
            sweep = 2.0 * math.pi if e[0] == "C" else abs(arc_angles(region, e)[3])
            n = max(n, int(round(MIN_CIRCLE_NODES * sweep / (2.0 * math.pi))))
        out.append(edge_samples(region, e, n))
    return np.concatenate(out)


def signed_area(P):
    Q = np.roll(P, -1, axis=0)
    return 0.5 * float(np.sum(P[:, 0] * Q[:, 1] - Q[:, 0] * P[:, 1]))


def convex_hull(P):
    """Andrew's monotone chain, counter-clockwise (Mesh3Dto2Dfacesgmsh's convexHull)."""
    pts = sorted(set(map(tuple, np.round(P, 15))))
    if len(pts) < 3:
        return np.array(pts)
    cross = lambda o, a, b: (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0.0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0.0:
            upper.pop()
        upper.append(p)
    return np.array(lower[:-1] + upper[:-1])


def pose(region):
    """Mesh3Dto2Dfacesgmsh's pose: the minimum-area rectangle along the axes, long side on
    x spanning [-1, 1], centred. Ties go to the shorter long side, as there."""
    hull = convex_hull(fine_samples(region, region.loops[0]))
    best = None
    for k in range(len(hull)):
        d = hull[(k + 1) % len(hull)] - hull[k]
        if not np.linalg.norm(d) > 0.0:
            continue
        d = d / np.linalg.norm(d)
        n = np.array([-d[1], d[0]])
        u, v = hull @ d, hull @ n
        w, ht = u.max() - u.min(), v.max() - v.min()
        area, long_side = w * ht, max(w, ht)
        if best is None or area < best[0] * (1 - 1e-9) or (area <= best[0] * (1 + 1e-9) and long_side < best[1]):
            centre = d * 0.5 * (u.min() + u.max()) + n * 0.5 * (v.min() + v.max())
            best = (area, long_side, d if w >= ht else n, centre)
    if best is None or not best[1] > 0.0:
        raise Reject("invalid", "a degenerate profile")
    _, long_side, X, centre = best
    Y = np.array([-X[1], X[0]])
    s = 2.0 / long_side
    return region.transformed(lambda p: s * np.array([np.dot(p - centre, X), np.dot(p - centre, Y)]), s)


# ----------------------------------------------------------------------------------------
# Checks, features, signature
# ----------------------------------------------------------------------------------------

def inside(p, P):
    """Point p strictly inside polygon P (even-odd rule)."""
    Q = np.roll(P, -1, axis=0)
    cond = (P[:, 1] > p[1]) != (Q[:, 1] > p[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        x = P[:, 0] + (p[1] - P[:, 1]) * (Q[:, 0] - P[:, 0]) / (Q[:, 1] - P[:, 1])
    return bool(np.count_nonzero(cond & (p[0] < x)) % 2)


def crossing(polys):
    """True when two segments of the loops' polylines meet anywhere but at the point two
    neighbours on one loop share. Sweep-and-prune on x, then the pairs all at once."""
    A, B, loop, idx, size = [], [], [], [], []
    for l, P in enumerate(polys):
        A.append(P)
        B.append(np.roll(P, -1, axis=0))
        loop.append(np.full(len(P), l))
        idx.append(np.arange(len(P)))
        size.append(np.full(len(P), len(P)))
    A, B = np.concatenate(A), np.concatenate(B)
    loop, idx, size = np.concatenate(loop), np.concatenate(idx), np.concatenate(size)
    lo, hi = np.minimum(A, B), np.maximum(A, B)
    order = np.argsort(lo[:, 0], kind="stable")
    xs = lo[order, 0]
    end = np.searchsorted(xs, hi[order, 0], side="right")
    count = np.maximum(end - np.arange(len(order)) - 1, 0)
    if count.sum() == 0:
        return False
    first = np.repeat(np.arange(len(order)), count)
    offset = np.arange(count.sum()) - np.repeat(np.cumsum(count) - count, count)
    i, j = order[first], order[first + 1 + offset]
    keep = (lo[i, 1] <= hi[j, 1]) & (lo[j, 1] <= hi[i, 1])
    same = loop[i] == loop[j]
    gap = np.abs(idx[i] - idx[j])
    keep &= ~(same & ((gap == 1) | (gap == size[i] - 1)))
    i, j = i[keep], j[keep]
    if len(i) == 0:
        return False
    p, q, r, s = A[i], B[i], A[j], B[j]
    cr = lambda o, a, b: (a[:, 0] - o[:, 0]) * (b[:, 1] - o[:, 1]) - (a[:, 1] - o[:, 1]) * (b[:, 0] - o[:, 0])
    d1, d2, d3, d4 = cr(r, s, p), cr(r, s, q), cr(p, q, r), cr(p, q, s)
    hit = (d1 * d2 <= 0.0) & (d3 * d4 <= 0.0)
    # All four zero: collinear, and the boxes overlapping (checked above) is meeting.
    return bool(hit.any())


def turns(P):
    """Signed turning at each vertex of closed polyline P, + to the left (BoundaryFeatures')."""
    u = P - np.roll(P, 1, axis=0)                 # into each vertex
    w = np.roll(P, -1, axis=0) - P                # out of it
    return np.arctan2(u[:, 0] * w[:, 1] - u[:, 1] * w[:, 0], np.sum(u * w, axis=1))


def features(polys):
    """BoundaryFeatures on the loops' polylines (region on the left): area, perimeter,
    corner angles, runs (length, turning), and the principal second moments of area."""
    area = sum(signed_area(P) for P in polys)
    perimeter, angles, runs = 0.0, [], []
    for P in polys:
        turn = turns(P)
        edge = np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1)
        perimeter += edge.sum()
        corner = np.abs(turn) > CORNER
        at = np.flatnonzero(corner)
        angles += list(np.degrees(math.pi - turn[at]))
        if len(at) <= 1:
            runs.append((edge.sum(), float(turn[~corner].sum())))
            continue
        for q in range(len(at)):
            a, b = at[q], at[(q + 1) % len(at)]
            span = np.arange(a, b) if b > a else np.concatenate([np.arange(a, len(P)), np.arange(0, b)])
            runs.append((float(edge[span].sum()), float(turn[span[1:]].sum())))
    # Second moments about the centroid (Green's theorem on each loop; holes are clockwise,
    # so they subtract).
    sx = sy = sxx = syy = sxy = 0.0
    for P in polys:
        x0, y0 = P[:, 0], P[:, 1]
        x1, y1 = np.roll(x0, -1), np.roll(y0, -1)
        c = x0 * y1 - x1 * y0
        sx += float(np.sum((x0 + x1) * c)) / 6.0
        sy += float(np.sum((y0 + y1) * c)) / 6.0
        sxx += float(np.sum((x0 * x0 + x0 * x1 + x1 * x1) * c)) / 12.0
        syy += float(np.sum((y0 * y0 + y0 * y1 + y1 * y1) * c)) / 12.0
        sxy += float(np.sum((x0 * y1 + 2 * x0 * y0 + 2 * x1 * y1 + x1 * y0) * c)) / 24.0
    cx, cy = sx / area, sy / area
    sxx, syy, sxy = sxx - area * cx * cx, syy - area * cy * cy, sxy - area * cx * cy
    mean, dev = 0.5 * (sxx + syy), math.sqrt(0.25 * (sxx - syy) ** 2 + sxy * sxy)
    return dict(area=area, perimeter=perimeter, angles=angles, runs=runs,
                moments=((mean - dev) / area ** 2, (mean + dev) / area ** 2))


def signature(f, loops):
    """Mesh3Dto2Dfacesgmsh's Signature: (loops, n angles, n runs) and the numbers sameShape
    compares, as one vector: perimeter, moments, angles, runs, turnings."""
    root = math.sqrt(f["area"])
    angles = sorted(f["angles"])
    runs = sorted(r[0] / root for r in f["runs"])
    turnings = sorted(math.degrees(r[1]) for r in f["runs"])
    key = (loops, len(angles), len(runs))
    return key, np.array([f["perimeter"] / root, *f["moments"], *angles, *runs, *turnings])


def same_shape(a, b, n_angles, n_runs):
    """Mesh3Dto2Dfacesgmsh's sameShape on two signature vectors with one key."""
    def near(x, y, rel, ab):
        return np.all(np.abs(x - y) <= ab + rel * np.maximum(np.abs(x), np.abs(y)))
    i = 3 + n_angles
    return (near(a[0], b[0], 3e-3, 0.0) and near(a[1:3], b[1:3], 5e-3, 0.0)
            and near(a[3:i], b[3:i], 0.0, 1.0) and near(a[i:i + n_runs], b[i:i + n_runs], 5e-3, 1e-3)
            and near(a[i + n_runs:], b[i + n_runs:], 0.0, 1.5))


# ----------------------------------------------------------------------------------------
# The .geo
# ----------------------------------------------------------------------------------------

def clean(s):
    """A name from the JSON, safe inside a // comment."""
    return "".join(ch if ch.isprintable() else " " for ch in str(s))[:80]


def geo_text(region, head):
    """The posed region as an OpenCASCADE .geo (see the docstring). Every Point comes before
    the first curve: a whole circle Circle(n) = {x, y, 0, r} makes its own point under the
    next free tag, which a Point written after it would then collide with."""
    used = sorted({e[k] for L in region.loops for e in L if e[0] != "C" for k in (1, 2)})
    tag = {v: n + 1 for n, v in enumerate(used)}
    f = lambda x: repr(float(x))
    lines = head + ['SetFactory("OpenCASCADE");', ""]
    for v in used:
        p = region.points[v]
        lines.append("Point(%d) = {%s, %s, 0};" % (tag[v], f(p[0]), f(p[1])))
    centres, curve, loop_tags = [], 0, []
    body = []
    for L in region.loops:
        tags = []
        for e in L:
            curve += 1
            if e[0] == "L":
                body.append("Line(%d) = {%d, %d};" % (curve, tag[e[1]], tag[e[2]]))
            elif e[0] == "A":
                c = len(used) + len(centres) + 1
                lines.append("Point(%d) = {%s, %s, 0};   // centre of curve %d" % (c, f(e[3][0]), f(e[3][1]), curve))
                body.append("Circle(%d) = {%d, %d, %d};" % (curve, tag[e[1]], c, tag[e[2]]))
                centres.append(c)
            else:
                body.append("Circle(%d) = {%s, %s, 0, %s};" % (curve, f(e[1][0]), f(e[1][1]), f(e[2])))
            tags.append(curve)
        loop_tags.append(tags)
    lines += [""] + body + [""]
    for k, tags in enumerate(loop_tags):
        lines.append("Curve Loop(%d) = {%s};" % (k + 1, ", ".join(map(str, tags))))
    lines.append("Plane Surface(1) = {%s};" % ", ".join(str(k + 1) for k in range(len(loop_tags))))
    if centres:
        lines += ["", "// The arc centres last: they are not part of the face (see the header).",
                  "Delete { Point{%s}; }" % ", ".join(map(str, centres))]
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------------------
# One model
# ----------------------------------------------------------------------------------------

def profiles_of(model):
    """(k, sketch name, profile id, extruded, loops) of every profile, sketches in sequence
    order, each sketch's profiles in the file's order."""
    ents = model["entities"]
    extruded = set()
    for e in ents.values():
        if e.get("type") == "ExtrudeFeature":
            for p in e.get("profiles", []):
                extruded.add((p["sketch"], p["profile"]))
    k = 0
    for item in model.get("sequence", []):
        if item.get("type") != "Sketch":
            continue
        sk = ents[item["entity"]]
        for pid, prof in sk.get("profiles", {}).items():
            yield k, sk.get("name", item["entity"]), pid, (item["entity"], pid) in extruded, prof["loops"]
            k += 1


def convert(raw_loops, h, min_run, min_width):
    """One profile: the posed region, its features and signature, or Reject."""
    region = pose(build_region(read_curves(raw_loops)))
    polys = [mesh_samples(region, L, h) for L in region.loops]
    if any(len(P) < 3 for P in polys):
        raise Reject("invalid", "a loop of fewer than three mesh edges")
    for P in polys[1:]:
        if not inside(P[0], polys[0]):
            raise Reject("invalid", "a hole outside the outer loop")
        if any(inside(P[0], Q) for Q in polys[1:] if Q is not P):
            raise Reject("invalid", "a hole inside another")
    if crossing(polys):
        raise Reject("invalid", "the boundary crosses itself at the mesh size")
    f = features(polys)
    if not f["area"] > 0.0:
        raise Reject("invalid", "no area")
    run = min(r[0] for r in f["runs"])
    width = 2.0 * f["area"] / f["perimeter"]
    if any(np.any(np.abs(turns(P)) > CUSP) for P in polys):
        raise Reject("not nice", "a cusp (a corner sharper than 1 degree)")
    if run < min_run * h:
        raise Reject("not nice", "a feature under %g mesh edges long" % min_run)
    if width < min_width * h:
        raise Reject("not nice", "under %g mesh edges thick" % min_width)
    return region, f, run / h, width / h


def model_task(args):
    """Everything one model's profiles give: a list of dicts, in profile order."""
    path, model_id, split, h, min_run, min_width, all_profiles = args
    try:
        with open(path) as fp:
            model = json.load(fp)
        found = list(profiles_of(model))
    except Exception as exc:
        return [dict(k=-1, status="invalid", why="unreadable JSON (%s)" % type(exc).__name__)]
    out = []
    for k, sketch, pid, extruded, loops in found:
        if not (extruded or all_profiles):
            continue
        rec = dict(k=k, sketch=clean(sketch), profile=clean(pid), extruded=extruded)
        try:
            region, f, run, width = convert(loops, h, min_run, min_width)
        except Reject as r:
            rec.update(status=r.kind, why=r.why)
            out.append(rec)
            continue
        except Exception as exc:   # a malformed profile is one rejection, not the run's end
            rec.update(status="invalid", why="unreadable profile (%s)" % type(exc).__name__)
            out.append(rec)
            continue
        kinds = Counter(e[0] for L in region.loops for e in L)
        key, sig = signature(f, len(region.loops))
        name = "%s_%d" % (model_id.split("/")[-1], k)
        head = ["%s%s: %s, profile %s (%s; DeepCAD split %s)." % (
                    HEADER, model_id, rec["sketch"], rec["profile"],
                    "extruded" if extruded else "not extruded", split),
                "// Written by py/deepcad_to_geo.py, which says how it is posed and why the arc centres go.",
                ""]
        rec.update(status="ok", name=name, key=key, sig=sig, geo=geo_text(region, head),
                   loops=len(region.loops), lines=kinds["L"], arcs=kinds["A"], circles=kinds["C"],
                   corners=len(f["angles"]), run_h=round(run, 2), width_h=round(width, 2))
        out.append(rec)
    return out


# ----------------------------------------------------------------------------------------
# The run
# ----------------------------------------------------------------------------------------

class Registry:
    """The shapes written so far, bucketed by (loops, n angles, n runs) and sorted by
    perimeter, so a new shape is compared only with those within sameShape's perimeter
    tolerance."""

    def __init__(self):
        self.buckets = {}

    def twin(self, key, sig):
        perims, sigs, names = self.buckets.get(key, ([], [], []))
        p = sig[0]
        lo, hi = bisect.bisect_left(perims, p * (1 - 3e-3)), bisect.bisect_right(perims, p / (1 - 3e-3))
        for n in range(lo, hi):
            if same_shape(sig, sigs[n], key[1], key[2]):
                return names[n]
        return None

    def add(self, key, sig, name):
        perims, sigs, names = self.buckets.setdefault(key, ([], [], []))
        n = bisect.bisect_left(perims, sig[0])
        perims.insert(n, sig[0])
        sigs.insert(n, sig)
        names.insert(n, name)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=DEFAULT_DATA, help="DeepCAD's data folder (cad_json/, train_val_test_split.json)")
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--np", type=int, default=100, help="make_meshes.sh's NP: the mesh size is h = 2/np")
    ap.add_argument("--min-run", type=float, default=3.0, help="shortest corner-to-corner run, in mesh edges")
    ap.add_argument("--min-width", type=float, default=4.0, help="thinnest 2A/P, in mesh edges")
    ap.add_argument("--all-profiles", action="store_true", help="also the sketch regions no feature extruded")
    ap.add_argument("--limit", type=int, default=0, help="only the first N models (by id)")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.data, "cad_json", "*", "*.json")))
    if not files:
        sys.exit("no DeepCAD JSON under %s/cad_json" % args.data)
    if args.limit:
        files = files[:args.limit]
    split = {}
    split_path = os.path.join(args.data, "train_val_test_split.json")
    if os.path.exists(split_path):
        with open(split_path) as fp:
            for name, ids in json.load(fp).items():
                for i in ids:
                    split[i] = name
    model_ids = [os.path.relpath(f, os.path.join(args.data, "cad_json"))[:-len(".json")] for f in files]

    os.makedirs(args.out, exist_ok=True)
    stale = 0
    for g in glob.glob(os.path.join(args.out, "*.geo")):
        with open(g) as fp:
            if fp.readline().startswith(HEADER):
                os.remove(g)
                stale += 1
    h = 2.0 / args.np
    print("%d models -> %s (h = 2/%d; %d old .geo removed)" % (len(files), args.out, args.np, stale), flush=True)

    reg, rows, copies = Registry(), [], Counter()
    status, why = Counter(), Counter()
    tasks = [(f, m, split.get(m, "none"), h, args.min_run, args.min_width, args.all_profiles)
             for f, m in zip(files, model_ids)]
    start = time.time()
    with Pool(args.jobs) as pool:
        for n, (model_id, recs) in enumerate(zip(model_ids, pool.imap(model_task, tasks, chunksize=16))):
            for rec in recs:
                if rec["status"] != "ok":
                    status[rec["status"]] += 1
                    why[(rec["status"], rec["why"])] += 1
                    continue
                twin = reg.twin(rec["key"], rec["sig"])
                if twin is not None:
                    status["duplicate"] += 1
                    copies[twin] += 1
                    continue
                reg.add(rec["key"], rec["sig"], rec["name"])
                copies[rec["name"]] += 1
                with open(os.path.join(args.out, rec["name"] + ".geo"), "w") as fp:
                    fp.write(rec["geo"])
                status["written"] += 1
                rows.append(dict(name=rec["name"], model=model_id, split=split.get(model_id, "none"),
                                 sketch=rec["sketch"], profile=rec["profile"], extruded=int(rec["extruded"]),
                                 **{k: rec[k] for k in ("loops", "lines", "arcs", "circles", "corners",
                                                        "run_h", "width_h")}))
            if (n + 1) % 5000 == 0 or n + 1 == len(tasks):
                print("  %d/%d models, %.0f s: %s" % (n + 1, len(tasks), time.time() - start,
                                                     ", ".join("%s %d" % kv for kv in sorted(status.items()))),
                      flush=True)

    with open(os.path.join(args.out, "manifest.csv"), "w", newline="") as fp:
        cols = ["name", "model", "split", "sketch", "profile", "extruded", "loops", "lines", "arcs",
                "circles", "corners", "run_h", "width_h", "copies"]
        w = csv.DictWriter(fp, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(dict(r, copies=copies[r["name"]]))
    total = sum(status.values())
    print("%d profiles: %s" % (total, ", ".join("%s %d" % kv for kv in sorted(status.items()))))
    for (kind, reason), count in why.most_common():
        print("  %-8s %6d  %s" % (kind, count, reason))
    print("%d .geo + manifest.csv in %s" % (status["written"], args.out))


if __name__ == "__main__":
    main()
