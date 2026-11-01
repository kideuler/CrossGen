"""
quarry.py -- QUARRY (src/QUARRY, CrossGen's C++ block extractor) on an extraction mesh's fields.

SQuadGen's last step turns the decoded fields into a quad layout (docs/squadgen_2d.md Sec. 4 step 8). The release does
it with QuadTools' QuadExtraction, which the 2026-10-29 round trip found fragile on decoded fields (spurious charts,
non-quad polygons) and, in 2D, unsound even on true ones (patch sides drawn as chords that leave the domain). QUARRY
reads the same four fields -- cdf, dcdf, offsetc (to the chart centre) and offset (to the subchart's corner) at every
triangle's centroid -- and returns CrossGen's own BlockDecomposition, which BlockQuadMesh meshes and TMOP smooths as
every other method's blocks are. Mesh.quarry() in the crossgen module is the call; this module wraps it for the
squadgen2d tools.

It needs the crossgen module and nothing from torch. That is deliberate, not tidiness: torch and crossgen each link
their own OpenMP runtime, and a process that has initialised torch's aborts when crossgen's TMOP starts its own. So the
round trip, which has torch loaded to decode, runs this file as a subprocess per extraction (run_quarry in
roundtrip.py), the way it runs QuadExtraction.

    python -m squadgen2d quarry DIR      # DIR/quarry_input.npz -> DIR/quarry_output.npz
    python py/squadgen/squadgen2d/quarry.py DIR

quarry_input.npz: V (n, 2), T (F, 3), cdf (F,), dcdf (F,), offset (F, 2), offsetc (F, 2), all in the mesh's units;
optional target (the mesh's target edge length; default 0.01 x the bounding-box diagonal, the CDF data's) and sweeps
(TMOP sweeps on the mesh; default 300).

quarry_output.npz: labels (F,) the chart of every triangle (-1 none); covers (the decomposition covers the domain:
every chart a quad, every side shared or on dS); vertices (n, 2), blocks (B, 4), edge_points (P, 2) and edge_start
(E + 1,) the macro edges as polylines; failed (charts that made no block); and, when the blocks cover, the mesh after
TMOP -- mesh_vertices, mesh_quads, inverted, min_sj -- and before it, inverted_raw and min_sj_raw; report (JSON).
"""

import json
import os
import sys
import time

import numpy as np

DEFAULT_SWEEPS = 300
DEFAULT_TARGET_REL = 0.01


def extract(V, T, cdf, dcdf, offset, offsetc, target=None, sweeps=DEFAULT_SWEEPS, **options):
    """QUARRY on per-triangle fields of the mesh (V, T). Returns the dict quarry_output.npz holds."""
    import crossgen

    t0 = time.time()
    m = crossgen.Mesh(np.asarray(V, float), np.asarray(T, np.int64))
    b = m.quarry(cdf=np.asarray(cdf, float).reshape(-1), dcdf=np.asarray(dcdf, float).reshape(-1),
                 offsetc=np.asarray(offsetc, float)[:, :2], offset=np.asarray(offset, float)[:, :2], **options)
    r = b.report
    out = {"labels": np.asarray(r["chart_labels"], np.int32), "covers": bool(b.valid),
           "vertices": np.asarray(b.vertices, float).reshape(-1, 2),
           "blocks": np.asarray(b.blocks, np.int64).reshape(-1, 4),
           "failed": int(r["failed_charts"]) + int(r["unmatched_sides"])}
    pts = [np.asarray(e, float).reshape(-1, 2) for e in b.edges]
    out["edge_points"] = np.concatenate(pts) if pts else np.zeros((0, 2))
    out["edge_start"] = np.cumsum([0] + [len(p) for p in pts])
    rep = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()
           if k not in ("chart_labels", "dual_labels")}
    if b.valid and len(out["blocks"]):
        if target is None:
            target = DEFAULT_TARGET_REL * float(np.linalg.norm(np.ptp(np.asarray(V, float), 0)))
        q = b.mesh(float(target))
        raw = q.quality
        out.update(inverted_raw=int(raw["inverted_quads"]), min_sj_raw=float(raw["min_scaled_jacobian"]))
        if sweeps:
            q.smooth(int(sweeps))
        qq = q.quality
        out.update(mesh_vertices=np.asarray(q.vertices, float), mesh_quads=np.asarray(q.quads, np.int64),
                   inverted=int(qq["inverted_quads"]), min_sj=float(qq["min_scaled_jacobian"]))
    rep["seconds_total"] = time.time() - t0
    out["report"] = json.dumps(rep)
    return out


def run_dir(d):
    with np.load(os.path.join(d, "quarry_input.npz")) as z:
        a = {k: z[k] for k in z.files}
    kw = {}
    if "target" in a:
        kw["target"] = float(a.pop("target"))
    if "sweeps" in a:
        kw["sweeps"] = int(a.pop("sweeps"))
    out = extract(a["V"], a["T"], a["cdf"], a["dcdf"], a["offset"], a["offsetc"], **kw)
    np.savez(os.path.join(d, "quarry_output.npz"), **out)
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if len(argv) != 1 or argv[0] in ("-h", "--help"):
        print("usage: python -m squadgen2d quarry DIR    (DIR/quarry_input.npz -> DIR/quarry_output.npz)")
        return 2
    out = run_dir(argv[0])
    rep = json.loads(out["report"])
    print(f"{len(out['blocks'])} blocks, covers {out['covers']}, "
          + (f"mesh {len(out['mesh_quads'])} quads, {out['inverted']} inverted, min scaled Jacobian "
             f"{out['min_sj']:.3f}" if "mesh_quads" in out else "no mesh")
          + f" ({rep['seconds_total']:.2f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
