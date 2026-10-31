"""
roundtrip.py -- the SQ-VAE round trip through QuadExtraction (docs/squadgen_2d.md Sec. 6 step 5).

The guide calls this the key VAE test: encode a layout, decode its fields on a triangulation of
Omega, run QuadExtraction on them, and get back the charts the layout had. Wall IoU and L1
(sq_vae.py's eval) only say how close the fields are; this says whether they are close in
the way extraction needs. A gap in a decoded wall is a small L1 error, but region growing leaks
through it and two charts become one (Sec. 5.1 step 2).

Every layout is extracted twice, from the same triangulation and with the same flags:

* from its true fields: our C++ CDF (BlockDecomposition.cdf(), the NPZ's split mesh) evaluated
  at the triangle centroids with Appendix B's formulas. This is the control, the Sec. 6 step 1
  check on our own data. Where it fails, the layout cannot be recovered at this resolution
  whatever SQ-VAE does: the usual cause is a sliver block narrower than `ringsize` rings
  (Sec. 7, "Extraction mesh fine enough");
* from its decoded fields: SQ-VAE's posterior mean decoded at the same centroids, as SQ-Diffuse's
  samples will be at inference (Sec. 4 steps 5-8).

The triangulation is the domain's own OBJ (the mesh every method was given), refined by midpoint
subdivision until its median edge is about h_rel x the bounding-box diagonal. It is independent
of the layout, so no mesh edge lies along a true chart side unless dOmega puts it there. The
release remeshes to 0.2 % of the diagonal (get_triangle_centers, Sec. 4 step 6), which is the
default here. Subdivision keeps every polygon corner as a vertex, and corners are what -a reads.

QuadExtraction's input is the release's (infer_sqdiffuse.py): a PLY of the triangulation at
z = 0 in OBJ coordinates, and extract_mesh.npz with float32 cdf (F, 1), dcdf (F, 1), offset and
offsetc (F, 3, the dual-chart and chart centre, z = 0) and invT (3, 4), face order as in the
PLY. Offsets are mapped back from the SQ-VAE frame by the scale only, as de_norm_pca does. The
flags default to the release's: --ringsize=8 -a 150 --div=0 -t 0.1 -s 3 -r 15. QuadExtraction
runs with --debug=true for one file: cdf_final_clustered_face.ply, the chart (cluster) of every
input triangle, as a colour per cluster. That is how the extracted charts are compared with the
true ones. Non-quad polygons left after edge collapse are counted from non_quad_polygons.ply,
and QuadQuality scores the output mesh.

Scoring, per layout and per field source. A true chart is matched when some extracted chart
covers it with area IoU > 0.5 (at most one can). The layout is recovered when every true chart
is matched, nothing else was extracted, and every polygon is a quad. dec_vs_gt matches the
decoded extraction against the true-field one instead, which separates SQ-VAE's errors from
what extraction loses anyway.

QuadTools must be built with the Appendix A patch (quadtools_2d.patch). On macOS, Apple clang
rejects the release's sources (an MSVC default argument in EdgeCollapseInfo, and omp.h even with
USE_OPENMP=OFF); Homebrew GCC builds them as they are:
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-16 -DCMAKE_CXX_COMPILER=g++-16
    cmake --build build -j8 --target QuadExtraction QuadQuality CDFGen

Run
---
    python -m squadgen2d roundtrip                             # val + 48 train layouts of best.pt
    python -m squadgen2d roundtrip --models M9_3 B1_1          # just these domains
    python -m squadgen2d roundtrip --split val --keep          # keep QuadExtraction's inputs and debug files
    python -m squadgen2d roundtrip --ckpt py/squadgen/runs/sq_vae_2d_oracle/best.pt --data py/squadgen/data/cdf_oracle

The split is the SQ-VAE's (the families its Geom-AE held out, recorded in the checkpoint), and the layouts come from the
manifest and their NPZs, not a field pack: extraction needs each layout's split mesh.

Outputs in --out (default: runs/sq_vae_2d/roundtrip):
* roundtrip.csv: one row per layout, both extractions;
* roundtrip_<split>.png: true charts, extraction from true fields, extraction from decoded fields;
* cases/<split>/<model>__<method>/: gt_quad.ply, dec_quad.ply (the extracted meshes) and
  labels.npz (true, gt and dec chart of every extraction-mesh triangle);
* log.txt.

Other scripts (SQ-Diffuse inference) can import extraction_mesh, write_extraction_input and
run_extraction.
"""

import argparse
import concurrent.futures as cf
import csv
import json
import os
import re
import shutil
import subprocess
import time

import numpy as np
import torch

from . import geometry as geo
from . import training as tr
from .domains import geom_ae_val_families
from .layout import METHODS, Layout, layout_fields, read_manifest
from .paths import DEFAULT_CDF, DEFAULT_SQ_VAE, REPO
from .sq_vae import decode_fields, encode_layout, load_sq_vae

DEFAULT_RUN = DEFAULT_SQ_VAE
DEFAULT_QUADTOOLS = os.environ.get(
    "QUADTOOLS", os.path.join(os.path.dirname(REPO), "SQUADGEN", "QuadTools4SquadGen", "build"))
# Mirrors infer_sqdiffuse.py's run_quad_postprocess (debug_script/extract_quad_from_cdf.sh).
RELEASE_FLAGS = {"ringsize": 8, "a": 150.0, "div": 0, "t": 0.1, "s": 3, "r": 15}
MATCH_IOU = 0.5


# =============================================================================
# 1. EXTRACTION MESH
# =============================================================================


def subdivide(V, T):
    """One 1-to-4 midpoint subdivision. Orientation and every original vertex are kept."""
    E = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    uniq, inv = np.unique(np.sort(E, 1), axis=0, return_inverse=True)
    mid = len(V) + inv.ravel().reshape(3, -1).T                        # (F, 3): mids of 01, 12, 20
    V = np.concatenate([V, 0.5 * (V[uniq[:, 0]] + V[uniq[:, 1]])])
    a, b, c = T.T
    ab, bc, ca = mid.T
    return V, np.concatenate([np.stack(t, 1) for t in ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])


def extraction_mesh(obj, h_rel):
    """The domain OBJ, counter-clockwise, without unused vertices, subdivided until its median
    edge is at most 1.5 h_rel x the bounding-box diagonal. Returns (V, T, levels)."""
    V, T = geo.read_obj(obj)
    flip = geo.signed_areas(V, T) < 0
    T[flip] = T[flip][:, ::-1]
    used, T = np.unique(T, return_inverse=True)
    V, T = V[used], T.reshape(-1, 3)
    diag = np.linalg.norm(np.ptp(V, 0))
    levels = 0
    while levels < 4 and np.median(np.linalg.norm(V[T[:, 1]] - V[T[:, 0]], axis=1)) > 1.5 * h_rel * diag:
        V, T = subdivide(V, T)
        levels += 1
    return V, T, levels


def face_adjacency(T):
    """(f, g) pairs of triangles sharing an edge, and the boundary edges (vertex pairs)."""
    E = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    f = np.tile(np.arange(len(T)), 3)
    key = np.sort(E, 1)
    order = np.lexsort((key[:, 1], key[:, 0]))
    key, f, E = key[order], f[order], E[order]
    same = (key[1:] == key[:-1]).all(1)
    pairs = np.stack([f[:-1][same], f[1:][same]], 1)
    lone = np.ones(len(key), bool)
    lone[:-1] &= ~same
    lone[1:] &= ~same
    return pairs, E[lone]


def locate(V, T, P, chunk=40000):
    """The triangle of (V, T) holding each point of P, and P's barycentric coordinates in it.
    A point just outside the mesh gets the candidate it is least outside of, with barycentrics
    that still reproduce it: the split mesh follows a curved dOmega by coarser chords than the
    extraction mesh does, and layout_fields clamps its subchart coordinates anyway. Bucket grid,
    each triangle in every cell its slightly grown bounding box touches. Returns (tri, bary,
    outside), outside = how far the point is outside its triangle in barycentric units."""
    A = V[T]
    lo, hi = A.min(1), A.max(1)
    cell = 2 * np.median((hi - lo).max(1))
    origin = V.min(0) - cell
    n = np.ceil((V.max(0) + cell - origin) / cell).astype(int) + 1
    i0 = np.floor((lo - 0.25 * cell - origin) / cell).astype(int)
    i1 = np.floor((hi + 0.25 * cell - origin) / cell).astype(int)
    cells, tris = [], []
    span = (i1 - i0).max(0) + 1
    for dx in range(span[0]):
        for dy in range(span[1]):
            m = (i0[:, 0] + dx <= i1[:, 0]) & (i0[:, 1] + dy <= i1[:, 1])
            cells.append((i0[m, 0] + dx) * n[1] + i0[m, 1] + dy)
            tris.append(np.nonzero(m)[0])
    cells, tris = np.concatenate(cells), np.concatenate(tris)
    order = np.argsort(cells, kind="stable")
    cells, tris = cells[order], tris[order]
    start = np.searchsorted(cells, np.arange(n[0] * n[1] + 1))
    tri = np.empty(len(P), np.int64)
    bary = np.empty((len(P), 3))
    outside = np.empty(len(P))
    for s in range(0, len(P), chunk):
        p = P[s:s + chunk]
        ij = np.clip(np.floor((p - origin) / cell).astype(int), 0, n - 1)
        c = ij[:, 0] * n[1] + ij[:, 1]
        cnt = start[c + 1] - start[c]
        k = max(1, cnt.max())
        cand = start[c][:, None] + np.arange(k)[None]
        valid = np.arange(k)[None] < cnt[:, None]
        cand = np.where(valid, tris[np.minimum(cand, len(tris) - 1)], 0)
        a, b, cc = A[cand, 0], A[cand, 1], A[cand, 2]
        den = geo.cross2(b - a, cc - a)
        l1 = geo.cross2(p[:, None] - a, cc - a) / den
        l2 = geo.cross2(b - a, p[:, None] - a) / den
        lam = np.stack([1 - l1 - l2, l1, l2], -1)
        score = np.where(valid & (den != 0), lam.min(-1), -np.inf)
        best = score.argmax(1)
        r = np.arange(len(p))
        tri[s:s + chunk], bary[s:s + chunk], outside[s:s + chunk] = cand[r, best], lam[r, best], -score[r, best]
    # A point in an empty cell (far outside, which should not happen) takes the nearest centroid.
    lost = ~np.isfinite(outside)
    if lost.any():
        cen = A.mean(1)
        for i in np.nonzero(lost)[0]:
            t = np.argmin(((cen - P[i]) ** 2).sum(1))
            a, b, cc = A[t]
            den = geo.cross2(b - a, cc - a)
            l1, l2 = geo.cross2(P[i] - a, cc - a) / den, geo.cross2(b - a, P[i] - a) / den
            tri[i], bary[i] = t, (1 - l1 - l2, l1, l2)
            outside[i] = -min(bary[i])
    return tri, bary, np.maximum(outside, 0)


# =============================================================================
# 2. QUADEXTRACTION I/O
# =============================================================================


def write_ply(path, V, faces, colours=None):
    """Binary little-endian PLY, z = 0. faces: (F, k) array of one polygon size."""
    faces = np.asarray(faces)
    head = ["ply", "format binary_little_endian 1.0", f"element vertex {len(V)}",
            "property float x", "property float y", "property float z",
            f"element face {len(faces)}", "property list uchar int vertex_indices"]
    fdt = [("n", "u1"), ("i", "<i4", (faces.shape[1],))]
    if colours is not None:
        head += ["property uchar red", "property uchar green", "property uchar blue"]
        fdt += [("c", "u1", (3,))]
    rec = np.empty(len(faces), fdt)
    rec["n"], rec["i"] = faces.shape[1], faces
    if colours is not None:
        rec["c"] = colours
    xyz = np.zeros((len(V), 3), "<f4")
    xyz[:, :V.shape[1]] = V
    with open(path, "wb") as f:
        f.write(("\n".join(head + ["end_header"]) + "\n").encode())
        f.write(xyz.tobytes())
        f.write(rec.tobytes())


PLY_TYPES = {"char": "i1", "int8": "i1", "uchar": "u1", "uint8": "u1", "short": "<i2", "int16": "<i2",
             "ushort": "<u2", "uint16": "<u2", "int": "<i4", "int32": "<i4", "uint": "<u4", "uint32": "<u4",
             "float": "<f4", "float32": "<f4", "double": "<f8", "float64": "<f8"}


def read_ply(path):
    """The vertices (V, 3), polygons (a list of index arrays) and face colours (F, 3) or None of
    a binary little-endian PLY as happly writes them (QuadTools' output and debug meshes)."""
    with open(path, "rb") as f:
        data = f.read()
    end = data.index(b"end_header\n") + len(b"end_header\n")
    elems = []
    for line in data[:end].decode().splitlines():
        w = line.split()
        if w[0] == "format" and w[1] != "binary_little_endian":
            raise ValueError(f"{path}: {w[1]} PLY not supported")
        if w[0] == "element":
            elems.append((w[1], int(w[2]), []))
        elif w[0] == "property":
            elems[-1][2].append(("list", PLY_TYPES[w[2]], PLY_TYPES[w[3]], w[4]) if w[1] == "list"
                                else ("scalar", PLY_TYPES[w[1]], w[2]))
    pos, out = end, {}
    for name, count, props in elems:
        lists = [p for p in props if p[0] == "list"]
        if not lists:
            rec = np.frombuffer(data, np.dtype([(p[2], p[1]) for p in props]), count, pos)
            pos += rec.nbytes
            out[name] = rec
            continue
        # One list (vertex_indices) first, then scalars, as happly writes faces. If every list
        # has the length of the first, the element is a fixed record and numpy reads it at once.
        (_, ctype, itype, lname), scalars = lists[0], [p for p in props if p[0] == "scalar"]
        k = int(np.frombuffer(data, ctype, 1, pos)[0]) if count else 0
        dt = np.dtype([("n", ctype), (lname, itype, (k,))] + [(p[2], p[1]) for p in scalars])
        rec = np.frombuffer(data, dt, count, pos) if pos + dt.itemsize * count <= len(data) else None
        if rec is not None and (rec["n"] == k).all():
            pos += rec.nbytes
            out[name] = (list(rec[lname]), rec)
            continue
        polys, cols, sdt = [], [], np.dtype([(p[2], p[1]) for p in scalars])
        for _ in range(count):
            k = int(np.frombuffer(data, ctype, 1, pos)[0])
            pos += np.dtype(ctype).itemsize
            polys.append(np.frombuffer(data, itype, k, pos).copy())
            pos += np.dtype(itype).itemsize * k
            cols.append(np.frombuffer(data, sdt, 1, pos)[0] if scalars else None)
            pos += sdt.itemsize
        out[name] = (polys, np.array(cols, sdt) if scalars else None)
    v = out["vertex"]
    V = np.stack([v["x"], v["y"], v["z"]], 1).astype(np.float64)
    polys, rec = out.get("face", ([], None))
    colours = None
    if rec is not None and rec.dtype.names and "red" in rec.dtype.names:
        colours = np.stack([rec["red"], rec["green"], rec["blue"]], 1)
    return V, polys, colours


def write_extraction_input(d, V, T, fields, centre, scale):
    """QuadExtraction's input as infer_sqdiffuse.py writes it: mesh.ply (z = 0, OBJ
    coordinates) and extract_mesh.npz. `fields` holds cdf, dcdf, offset and offsetc per triangle
    in the SQ-VAE frame; offsets go back by 1/scale (de_norm_pca's rotation-only map) and invT
    is the frame's inverse, x_obj = x / scale + centre."""
    os.makedirs(d, exist_ok=True)
    write_ply(os.path.join(d, "mesh.ply"), V, T)
    F = len(T)
    out = {"cdf": fields["cdf"].reshape(F, 1), "dcdf": fields["dcdf"].reshape(F, 1)}
    for k in ("offset", "offsetc"):
        out[k] = np.concatenate([fields[k][:, :2] / scale, np.zeros((F, 1))], 1)
    invT = np.zeros((3, 4))
    invT[:, :3] = np.eye(3) / scale
    invT[:2, 3] = centre
    out["invT"] = invT
    np.savez(os.path.join(d, "extract_mesh.npz"), **{k: v.astype(np.float32) for k, v in out.items()})


def run_extraction(bin_dir, d, flags, timeout=900):
    """QuadExtraction on d/mesh.ply + d/extract_mesh.npz -> d/quad.ply, then QuadQuality on it.
    Returns a dict: rc, seconds, the cluster (chart) of every input triangle, the number of
    non-quad polygons, the output mesh and QuadQuality's scores. Debug mode only adds files;
    the output is byte-identical (checked on M9_3 and B1_1). Non-quad polygons are read from the
    log: without --improve they are not fixed, and the -s subdivision splits an n-gon into n
    quads, so the output mesh is all quads either way and non_quad_polygons.ply (written after
    subdivision) is never there."""
    cmd = [os.path.join(bin_dir, "QuadExtraction"), "-i", os.path.join(d, "mesh.ply"),
           "-f", os.path.join(d, "extract_mesh.npz"), "-o", os.path.join(d, "quad.ply"),
           f"--ringsize={flags['ringsize']}", "-a", str(flags["a"]), f"--div={flags['div']}",
           "-t", str(flags["t"]), "-s", str(flags["s"]), "-r", str(flags["r"]), "--debug=true", "--verbose=true"]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        rc, log = p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        rc, log = -1, "timeout"
    res = {"rc": rc, "seconds": time.time() - t0, "labels": None, "nonquad": 0, "quad": None, "quality": {}}
    with open(os.path.join(d, "extraction.log"), "w") as f:
        f.write(log)
    lab = os.path.join(d, "cdf_final_clustered_face.ply")
    if os.path.exists(lab):
        _, _, col = read_ply(lab)
        res["labels"] = np.unique(col.astype(np.int64) @ np.array([65536, 256, 1]), return_inverse=True)[1].ravel()
    nq = re.search(r"Warning: (\d+) non-quad polygons exist", log)
    res["nonquad"] = int(nq.group(1)) if nq else 0
    q = os.path.join(d, "quad.ply")
    if rc == 0 and os.path.exists(q):
        V, polys, _ = read_ply(q)
        res["quad"] = (V[:, :2], polys)
        subprocess.run([os.path.join(bin_dir, "QuadQuality"), "-i", q, "-j", d], capture_output=True,
                       timeout=timeout)
        js = os.path.join(d, "quad_quadquality.json")
        if os.path.exists(js):
            with open(js) as f:
                res["quality"] = json.load(f)[0]
    return res


# =============================================================================
# 3. SCORING
# =============================================================================


def match(ref, found, area, n_ref=None):
    """Area IoU of every reference region with every found one. Returns matched (reference
    regions with a found region of IoU > MATCH_IOU, which is then unique), the area fraction
    they cover, the mean best IoU, and the number of found regions."""
    n_ref = ref.max() + 1 if n_ref is None else n_ref
    nf = found.max() + 1
    inter = np.bincount(ref * nf + found, weights=area, minlength=n_ref * nf).reshape(n_ref, nf)
    a_ref, a_found = inter.sum(1), inter.sum(0)
    iou = inter / np.maximum(a_ref[:, None] + a_found[None] - inter, 1e-300)
    best = iou.max(1)
    ok = best > MATCH_IOU
    return {"found": int((a_found > 0).sum()), "matched": int(ok.sum()),
            "area_matched": float(a_ref[ok].sum() / a_ref.sum()), "miou": float(best.mean())}


def mesh_validity(quad, V, T):
    """Whether the output mesh is a mesh of Omega: (quads inverted against the mesh's own
    orientation, area fraction of quads whose centroid is outside Omega). In 2D QuadExtraction
    can place a patch side as a straight chord between its corners where the true side follows
    a curved dOmega, and then the quads cross a hole or leave the domain."""
    QV, polys = quad
    sa = np.array([0.5 * np.sum(QV[p, 0] * QV[np.roll(p, -1), 1] - QV[np.roll(p, -1), 0] * QV[p, 1]) for p in polys])
    inverted = int((sa * np.sign(sa.sum()) <= 0).sum())
    cen = np.array([QV[p].mean(0) for p in polys])
    _, _, out = locate(V, T, cen)
    return inverted, float(np.abs(sa)[out > 0.05].sum() / np.abs(sa).sum())


def score(res, true, n_true, area, V, T):
    """`ok`: the charts were recovered (every true chart matched, nothing else found, every
    polygon a quad). `valid`: the output mesh has no inverted quad and stays inside Omega."""
    if res["labels"] is None or len(res["labels"]) != len(true):
        return {"rc": res["rc"], "found": 0, "matched": 0, "area_matched": 0.0, "miou": 0.0,
                "nonquad": res["nonquad"], "ok": False, "inverted": -1, "outside_area": float("nan"),
                "valid": False, "loop": float("nan"), "seconds": res["seconds"]}
    m = match(true, res["labels"], area, n_true)
    q = res["quality"]
    loop = min(q["Fratio"], q["Eratio"]) if q else float("nan")
    ok = res["rc"] == 0 and res["quad"] is not None and m["found"] == n_true == m["matched"] and res["nonquad"] == 0
    inv, outside = mesh_validity(res["quad"], V, T) if res["quad"] is not None else (-1, float("nan"))
    return dict(rc=res["rc"], **m, nonquad=res["nonquad"], ok=bool(ok), inverted=inv,
                outside_area=round(outside, 4), valid=bool(inv == 0 and outside < 1e-3), loop=loop,
                seconds=res["seconds"])


# =============================================================================
# 4. ROUND TRIP
# =============================================================================


def prepare(r, model, geom, cfg, dev, h_rel, work):
    """Both extraction inputs of one layout, and what scoring needs."""
    L = Layout(r["npz"])
    with np.load(r["npz"]) as z:
        chart_of_tri = z["subchart_block"][z["face2quad"].ravel()][np.arange(len(L.T)) // 2].astype(np.int64)
        n_true = int(z["subchart_block"].max()) + 1
    V, T, levels = extraction_mesh(r["path"], h_rel)
    P = (V[T].mean(1) - L.centre) * L.scale                            # centroids, SQ-VAE frame
    tri, bary, outside = locate(L.V, L.T, P)
    gt = layout_fields(L, L.V, tri, bary)
    z, _, _ = encode_layout(model, geom, cfg, r["npz"], dev)
    dec = decode_fields(model, z, P)
    for kind, fields in (("gt", gt), ("dec", dec)):
        write_extraction_input(os.path.join(work, kind), V, T, fields, L.centre, L.scale)
    return {"V": V, "T": T, "levels": levels, "true": chart_of_tri[tri], "n_true": n_true,
            "area": np.abs(geo.signed_areas(V, T)), "outside": float((outside > 1e-6).mean()), "work": work}


def finish(r, case, futs, out_dir, keep):
    res = {k: f.result() for k, f in futs.items()}
    row = {k: r[k] for k in ("model", "method", "split", "family", "path")}
    row.update(charts=case["n_true"], faces=len(case["T"]), levels=case["levels"], outside=round(case["outside"], 4))
    for k in ("gt", "dec"):
        row.update({f"{k}_{a}": v for a, v in score(res[k], case["true"], case["n_true"], case["area"], case["V"], case["T"]).items()})
    if res["gt"]["labels"] is not None and res["dec"]["labels"] is not None:
        m = match(res["gt"]["labels"], res["dec"]["labels"], case["area"])
        row["dec_vs_gt"] = round(m["matched"] / max(1, res["gt"]["labels"].max() + 1), 4)
    else:
        row["dec_vs_gt"] = 0.0
    d = os.path.join(out_dir, "cases", r["split"], f"{r['model']}__{r['method']}")
    os.makedirs(d, exist_ok=True)
    labels = {"true": case["true"].astype(np.int16)}
    for k in ("gt", "dec"):
        if res[k]["labels"] is not None:
            labels[k] = res[k]["labels"].astype(np.int16)
        q = os.path.join(case["work"], k, "quad.ply")
        if os.path.exists(q):
            shutil.copy(q, os.path.join(d, f"{k}_quad.ply"))
    np.savez_compressed(os.path.join(d, "labels.npz"), **labels)
    if keep:
        for k in ("gt", "dec"):
            dst = os.path.join(d, k)
            shutil.rmtree(dst, ignore_errors=True)
            shutil.move(os.path.join(case["work"], k), dst)
    shutil.rmtree(case["work"], ignore_errors=True)
    return row


def select_rows(args, val_fams):
    """The layouts to run: the SQ-VAE's val layouts and the 48 training layouts its train eval drew."""
    rows = read_manifest(args.data, args.methods)
    for r in rows:
        r["path"] = os.path.join(REPO, r["path"])
        if val_fams is not None:
            r["split"] = "val" if r["family"] in val_fams else "train"
    out = []
    if "val" in args.split:
        out += [r for r in rows if r["split"] == "val"]
    if "train" in args.split:
        # The 48 training layouts sq_vae.py's train eval used (same seed, same draw).
        train_rows = [r for r in rows if r["split"] == "train"]
        rng = np.random.default_rng(args.seed)
        out += [train_rows[i] for i in rng.choice(len(train_rows), min(args.train_n, len(train_rows)), replace=False)]
    if args.models:
        out = [r for r in out if r["model"] in args.models]
    return out


# =============================================================================
# 5. REPORT
# =============================================================================


def summarise(rows, log):
    def pct(a, b):
        return f"{a}/{b} ({100 * a / max(1, b):.0f}%)"

    for split in ("val", "train"):
        rs = [r for r in rows if r["split"] == split]
        if not rs:
            continue
        log(f"{split}: {len(rs)} layouts")
        for name, sel in (("all", lambda r: True), ("1 chart", lambda r: r["charts"] == 1),
                          ("2-9 charts", lambda r: 2 <= r["charts"] <= 9), ("10+ charts", lambda r: r["charts"] >= 10)):
            s = [r for r in rs if sel(r)]
            if not s:
                continue
            gt_ok = [r for r in s if r["gt_ok"]]
            log(f"  {name:10s} n={len(s):3d}  recovered: true fields {pct(len(gt_ok), len(s))}, "
                f"decoded {pct(sum(r['dec_ok'] for r in s), len(s))}, decoded where true fields recover "
                f"{pct(sum(r['dec_ok'] for r in gt_ok), len(gt_ok))}  |  charts matched: true "
                f"{sum(r['gt_matched'] for r in s)}/{sum(r['charts'] for r in s)}, decoded "
                f"{sum(r['dec_matched'] for r in s)}  |  median dec_vs_gt {np.median([r['dec_vs_gt'] for r in s]):.2f}  |  "
                f"valid mesh: true {sum(r['gt_valid'] for r in s)}, decoded {sum(r['dec_valid'] for r in s)}")
        fails = [r for r in rs if r["gt_ok"] and not r["dec_ok"]]
        if fails:
            log("  decoded fails where true fields recover: " + ", ".join(
                f"{r['model']}/{r['method']} ({r['dec_matched']}/{r['charts']} matched, {r['dec_found']} found"
                + (f", {r['dec_nonquad']} non-quad" if r["dec_nonquad"] else "") + ")" for r in fails))


# Light tints of the dataviz skill's categorical slots 1-6: region fills, where only adjacent
# regions must differ, so six colours (smallest-last greedy colouring of a planar map) suffice.
REGION_FILLS = ("#a9c9ef", "#f5b39a", "#97d9c1", "#f6d58f", "#f3c4d6", "#8fcb8f")


def region_colours(labels, pairs):
    """A fill index per region such that regions sharing an edge differ: greedy colouring in
    smallest-last order, which needs at most six colours on a planar adjacency graph."""
    n = labels.max() + 1
    a, b = labels[pairs[:, 0]], labels[pairs[:, 1]]
    keep = a != b
    adj = [set() for _ in range(n)]
    for u, v in zip(a[keep], b[keep]):
        adj[u].add(v)
        adj[v].add(u)
    deg = {u: len(adj[u]) for u in range(n)}
    order, left = [], set(range(n))
    while left:
        u = min(left, key=lambda x: deg[x])
        order.append(u)
        left.remove(u)
        for v in adj[u]:
            if v in left:
                deg[v] -= 1
    col = np.zeros(n, int)
    done = set()
    for u in reversed(order):
        used = {col[v] for v in adj[u] if v in done}
        col[u] = next(c for c in range(len(REGION_FILLS) + n) if c not in used)
        done.add(u)
    return col % len(REGION_FILLS)


def plot_cases(rows, args, out, title, k=8):
    """A row per layout: the true charts, the extraction from the true fields and the one from
    the decoded fields, each filled by region. The extracted output quad mesh (-s 3 -r 15) is
    drawn over its regions, the true chart sides over the true ones. Rows: the decoded failures
    first (lowest decoded IoU, one per domain), then the recovered layouts with the most charts."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap

    seen, picked = set(), []
    fails = sorted((r for r in rows if not r["dec_ok"]), key=lambda r: (not r["gt_ok"], r["dec_miou"]))
    good = sorted((r for r in rows if r["dec_ok"]), key=lambda r: -r["charts"])
    for pool, n in ((fails, k // 2), (good, k)):
        for r in pool:
            if len(picked) < n and r["model"] not in seen:
                seen.add(r["model"])
                picked.append(r)
    if not picked:
        return
    cmap = ListedColormap(REGION_FILLS)
    meshes = [extraction_mesh(r["path"], args.h_rel) for r in picked]
    aspect = [float(np.clip(np.ptp(V[:, 1]) / np.ptp(V[:, 0]), 0.3, 1.3)) for V, _, _ in meshes]
    fig, axes = plt.subplots(len(picked), 3, figsize=(9.6, sum(3.0 * a + 0.5 for a in aspect) + 0.6), squeeze=False,
                             gridspec_kw={"height_ratios": [a + 0.17 for a in aspect]})
    for row, r, (V, T, _) in zip(axes, picked, meshes):
        d = os.path.join(out, "cases", r["split"], f"{r['model']}__{r['method']}")
        lab = np.load(os.path.join(d, "labels.npz"))
        pairs, bnd = face_adjacency(T)
        heads = {"true": f"{r['charts']} true chart" + ("" if r["charts"] == 1 else "s")}
        for k_ in ("gt", "dec"):
            heads[k_] = (f"{r[k_ + '_found']} found, {r[k_ + '_matched']}/{r['charts']} matched"
                         + (f", {r[k_ + '_nonquad']} non-quad" if r[k_ + "_nonquad"] else "")
                         + ("  [recovered]" if r[k_ + "_ok"] else "")
                         + ("" if r[k_ + "_valid"] else "\nmesh invalid: " + (f"{r[k_ + '_inverted']} inverted, " if r[k_ + "_inverted"] else "")
                            + f"{100 * r[k_ + '_outside_area']:.1f}% of area outside"))
        for ax, k_ in zip(row, ("true", "gt", "dec")):
            if k_ in lab.files:
                L_ = lab[k_].astype(np.int64)
                ax.tripcolor(V[:, 0], V[:, 1], T, facecolors=region_colours(L_, pairs).astype(float)[L_], cmap=cmap,
                             vmin=-0.5, vmax=len(REGION_FILLS) - 0.5, shading="flat", edgecolors="none", rasterized=True)
                if k_ == "true":
                    side = L_[pairs[:, 0]] != L_[pairs[:, 1]]
                    E = shared_edges(T, pairs[side])
                    segs = np.concatenate([V[E], V[bnd]])
                    ax.add_collection(LineCollection(segs, colors="#1a1a19", linewidths=0.7))
            q = os.path.join(d, f"{k_}_quad.ply")
            if k_ != "true" and os.path.exists(q):
                QV, polys, _ = read_ply(q)
                segs = [QV[[p[i], p[(i + 1) % len(p)]], :2] for p in polys for i in range(len(p))]
                ax.add_collection(LineCollection(segs, colors="#1a1a19", linewidths=0.25, alpha=0.7))
            ax.set_title(heads[k_], fontsize=7.5, color="#52514e")
            ax.set_xlim(V[:, 0].min(), V[:, 0].max())
            ax.set_ylim(V[:, 1].min(), V[:, 1].max())
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
        row[0].set_ylabel(f"{r['model']} {r['method']}", fontsize=8, color="#1a1a19")
    for ax, c in zip(axes[0], ("true layout", "QuadExtraction of true fields", "QuadExtraction of decoded fields")):
        ax.set_title(c + "\n" + ax.get_title(), fontsize=8, color="#1a1a19")
    fig.suptitle(title, y=0.997, fontsize=9, color="#1a1a19")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(os.path.join(out, f"roundtrip_{rows[0]['split']}.png"), dpi=120, facecolor="#fcfcfb")
    plt.close(fig)


def shared_edges(T, pairs):
    """The vertex pair of the edge each (f, g) pair of adjacent triangles shares."""
    a, b = T[pairs[:, 0]], T[pairs[:, 1]]
    return a[(a[:, :, None] == b[:, None, :]).any(2)].reshape(-1, 2)


# =============================================================================
# 6. MAIN
# =============================================================================


def main(argv=None):
    ap = argparse.ArgumentParser(prog="squadgen2d roundtrip", description=__doc__.split("\n")[1].strip())
    ap.add_argument("--ckpt", default=os.path.join(DEFAULT_RUN, "best.pt"), help="sq_vae.py checkpoint")
    ap.add_argument("--data", default=DEFAULT_CDF, help="a cdf-data output directory (manifest.csv + npz/)")
    ap.add_argument("--out", default=os.path.join(DEFAULT_RUN, "roundtrip"))
    ap.add_argument("--quadtools", default=DEFAULT_QUADTOOLS, help="QuadTools build dir (QuadExtraction, QuadQuality)")
    ap.add_argument("--split", nargs="+", default=["val", "train"], choices=["val", "train"])
    ap.add_argument("--train_n", type=int, default=48, help="training layouts: sq_vae.py's train-eval draw")
    ap.add_argument("--methods", nargs="+", default=list(METHODS), choices=METHODS)
    ap.add_argument("--models", nargs="+", help="only these domains")
    ap.add_argument("--h_rel", type=float, default=0.002, help="extraction-mesh edge / bbox diagonal (release: 0.002)")
    ap.add_argument("--ringsize", type=int, default=RELEASE_FLAGS["ringsize"])
    ap.add_argument("--angle", type=float, default=RELEASE_FLAGS["a"], help="QuadExtraction -a")
    ap.add_argument("--div", type=int, default=RELEASE_FLAGS["div"])
    ap.add_argument("--threshold", type=float, default=RELEASE_FLAGS["t"], help="QuadExtraction -t")
    ap.add_argument("--subdiv", type=int, default=RELEASE_FLAGS["s"], help="QuadExtraction -s")
    ap.add_argument("--smooth", type=int, default=RELEASE_FLAGS["r"], help="QuadExtraction -r")
    ap.add_argument("--jobs", type=int, default=4, help="QuadExtraction processes at once")
    ap.add_argument("--keep", action="store_true", help="keep QuadExtraction's inputs and debug files per case")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", default="auto", help="cuda, cuda:N, mps or cpu (default: the first there is)")
    args = ap.parse_args(argv)

    for b in ("QuadExtraction", "QuadQuality"):
        if not os.access(os.path.join(args.quadtools, b), os.X_OK):
            raise SystemExit(f"no {b} in {args.quadtools}: build QuadTools (see this file's docstring) or pass --quadtools")
    flags = {"ringsize": args.ringsize, "a": args.angle, "div": args.div, "t": args.threshold,
             "s": args.subdiv, "r": args.smooth}
    log = tr.RunLog(args.out)
    dev = tr.pick_device(args.device)
    model, geom, cfg = load_sq_vae(args.ckpt, dev)
    ck = torch.load(args.ckpt, map_location="cpu", weights_only=False)
    step = ck["step"]
    if cfg.get("val_families") is not None:
        val_fams, how = set(cfg["val_families"]), "recorded in the SQ-VAE checkpoint"
    else:
        val_fams, how = geom_ae_val_families(ck["geom_config"])
    rows = select_rows(args, val_fams)
    log(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  {args.ckpt} (step {step})  {len(rows)} layouts  "
        f"h_rel={args.h_rel}  flags " + " ".join(f"{k}={v}" for k, v in flags.items()))
    log(f"split: the Geom-AE's held-out families, {how}" if val_fams is not None
        else f"split: the manifest's split column (the Geom-AE's held-out families: {how})")
    work_root = os.path.join(args.out, "_work")
    os.makedirs(work_root, exist_ok=True)
    pool = cf.ThreadPoolExecutor(args.jobs)
    pending, results, t0 = [], [], time.time()

    def drain(n):
        while len(pending) > n:
            r, case, futs = pending.pop(0)
            row = finish(r, case, futs, args.out, args.keep)
            results.append(row)
            log(f"[{len(results):3d}/{len(rows)}] {row['split']:5s} {row['model']:10s} {row['method']:8s} "
                f"charts {row['charts']:3d}  faces {row['faces']:6d}  "
                f"true: {row['gt_found']:3d} found {row['gt_matched']:3d} matched {row['gt_nonquad']} nq "
                f"{'OK ' if row['gt_ok'] else '-- '}{'' if row['gt_valid'] else '(invalid mesh) '} decoded: {row['dec_found']:3d} found {row['dec_matched']:3d} "
                f"matched {row['dec_nonquad']} nq {'OK ' if row['dec_ok'] else '-- '}{'' if row['dec_valid'] else '(invalid mesh) '} dec_vs_gt {row['dec_vs_gt']:.2f}  "
                f"({time.time() - t0:.0f} s)")

    for i, r in enumerate(rows):
        work = os.path.join(work_root, f"{i:04d}")
        case = prepare(r, model, geom, cfg, dev, args.h_rel, work)
        futs = {k: pool.submit(run_extraction, args.quadtools, os.path.join(work, k), flags) for k in ("gt", "dec")}
        pending.append((r, case, futs))
        drain(2 * args.jobs)
    drain(0)
    pool.shutdown()
    shutil.rmtree(work_root, ignore_errors=True)

    with open(os.path.join(args.out, "roundtrip.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0]))
        w.writeheader()
        w.writerows(results)
    summarise(results, log)
    for split in args.split:
        rs = [r for r in results if r["split"] == split]
        if rs:
            plot_cases(rs, args, args.out, f"SQ-VAE round trip through QuadExtraction, {split} layouts, "
                                            f"best.pt (step {step})")
    log(f"done in {(time.time() - t0) / 60:.1f} min -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
