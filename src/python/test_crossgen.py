"""End-to-end check of the crossgen module (ctest Python_Module).

    PYTHONPATH=build/python python3 src/python/test_crossgen.py data/meshes/singlemat/geom001.obj

On the one model given, every Mesh method is taken through the whole chain --
a block decomposition, a quad mesh on it at h = 0.05, a few TMOP sweeps -- and
what is asserted is structural: the blocks exist and index what they claim to,
the mesh has elements and no fold after smoothing, and the arrays have the
shapes and types the docstrings promise, every layout is valid and its three
decomposition qualities are in range. The Shape-DNA is checked against
itself (ascending, positive, and the same through a Mesh rebuilt from the
arrays), the boundary features against what a quarter disk owes and against a
mirrored copy, and the keyword errors are checked for being errors. How good a
layout is, is not asserted here: that is the drivers' business, and their
ctests'.
"""
import sys
import time

import numpy as np

import crossgen


def check(cond, what):
    if not cond:
        raise AssertionError(what)
    print(f"  ok  {what}")


def expect(exc, fn, what):
    try:
        fn()
    except exc as e:
        print(f"  ok  {what}: {type(e).__name__}: {e}")
        return
    raise AssertionError(f"{what}: no {exc.__name__} raised")


def main(path):
    print(f"crossgen {crossgen.__version__}, methods {crossgen.methods}")

    # ---- Mesh ---------------------------------------------------------------
    m = crossgen.load(path)
    print(m)
    V, T, M = m.vertices, m.triangles, m.materials
    check(V.dtype == np.float64 and V.shape == (m.num_vertices, 2), "vertices are (n, 2) float64")
    check(T.dtype == np.int64 and T.shape == (m.num_triangles, 3), "triangles are (m, 3) int64")
    check(M.shape == (m.num_triangles,), "one material per triangle")
    check(T.min() >= 0 and T.max() < m.num_vertices, "triangle indices are vertices")
    check(m.boundary_edges.shape[1] == 2 and len(m.boundary_edges) > 0, "boundary edges")

    # Rebuilt from its own arrays, with every triangle turned clockwise: the
    # constructor turns them back, so the spectrum cannot tell the two apart.
    m2 = crossgen.Mesh(V, T[:, ::-1], M)
    check(m2.num_triangles == m.num_triangles, "Mesh(vertices, triangles, materials) round trip")
    expect(ValueError, lambda: crossgen.Mesh(V, T + m.num_vertices), "an index out of range")
    expect(FileNotFoundError, lambda: crossgen.load(path + ".missing"), "a missing file")

    # ---- Shape-DNA ------------------------------------------------------------
    e = m.shape_dna(count=12, degree=2)
    check(e.shape == (12,) and np.all(e > 0) and np.all(np.diff(e) >= 0),
          f"shape_dna: 12 positive ascending values, first {e[0]:.4f}")
    e2, rep = m2.shape_dna(count=12, degree=2, report=True)
    check(np.allclose(e, e2, rtol=1e-9), "the same DNA from the rebuilt mesh")
    check(rep["converged"] and len(rep["eigenvalues"]) >= 12, "shape_dna report")
    w = m.shape_dna(count=12, degree=2, normalization="weyl_ratio")
    check(np.allclose(w, e / (4 * np.pi * np.arange(1, 13)), rtol=1e-12),
          f"weyl_ratio is the area-normalised DNA over 4 pi k: {w[0]:.4f} ... {w[-1]:.4f}")

    # ---- boundary features ----------------------------------------------------
    # geom001 is a quarter disk: three right angles, no hole, and T3 leaves any
    # quad layout of it one quarter turn short.
    f, rep = m.boundary_features(report=True)
    check(f["corners"] == 3 and f["corners_one_block"] == 3 and f["holes"] == 0,
          f"boundary_features: {f['corners']} right angles, {f['holes']} holes")
    check(abs(f["minimum_defect"] - 1.0) < 0.05, f"minimum_defect {f['minimum_defect']:.3f}: one quarter turn")
    check(rep["corners"].shape == (3, 2) and len(rep["loops"]) == 1, "boundary_features report")
    mirrored = crossgen.Mesh(V * np.array([-1.0, 1.0]), T, M).boundary_features()
    check(all(np.isclose(mirrored[k], f[k], rtol=1e-9, atol=1e-12) for k in f),
          "the same features from a mirrored copy")
    expect(TypeError, lambda: m.boundary_features(corner_angel=10.0), "a misspelt boundary_features keyword")

    # ---- options and their errors ---------------------------------------------
    for name in crossgen.methods:
        for stage in ("method", "mesh", "smooth"):
            check(isinstance(crossgen.options(name, stage), dict), f"options('{name}', '{stage}')")
    check(crossgen.options("meridian")["dual_mbo_gamma"] == 10.0, "a default read through")
    check(crossgen.options("zipline", "smooth")["quadrature"] == "corners", "an enum reads by name")
    expect(TypeError, lambda: m.zipline(feild_seed=3), "a misspelt keyword")
    expect(TypeError, lambda: m.shape_dna(count=2.5), "a float for an int")
    expect(ValueError, lambda: m.shape_dna(boundary="robin"), "an unknown enum name")
    expect(ValueError, lambda: crossgen.options("zipline", "tmop"), "an unknown stage")

    # ---- the five methods -------------------------------------------------------
    for name in crossgen.methods:
        t0 = time.time()
        b = getattr(m, name)()
        print(f"{b}  [{time.time() - t0:.1f} s]")
        nb = b.num_blocks
        check(nb > 0, f"{name}: {nb} block(s)")
        check(b.blocks.shape == (nb, 4) and b.block_edges.shape == (nb, 4), f"{name}: blocks are (b, 4)")
        check(b.blocks.max() < len(b.vertices), f"{name}: corners index macrovertices")
        check(b.block_edges.max() < len(b.edges), f"{name}: sides index macro edges")
        check(0.0 < b.coverage <= 1.0 + 1e-9, f"{name}: coverage {b.coverage:.3f}")
        check(isinstance(b.report, dict) and isinstance(b.options, dict), f"{name}: report and options")
        # Qualities in (0, 1], higher better, on a valid layout (NaN on one
        # that is not). regularity is measured against the quarter turn the
        # quarter disk owes, so a layout that owes no more reaches 1.
        check(b.valid, f"{name}: a valid decomposition")
        reg, ang, cho = b.regularity(), b.angle_quality(), b.chord_quality()
        check(0.0 < reg <= 1.0 and 0.0 < ang <= 1.0 and 0.0 < cho <= 1.0,
              f"{name}: regularity {reg:.3f}, angle quality {ang:.3f}, chord quality {cho:.3f}")

        q = b.mesh(0.05)
        print(f"  {q}")
        check(q.num_quads > 0 and q.quads.shape == (q.num_quads, 4), f"{name}: {q.num_quads} quads")
        check(q.scaled_jacobian.shape == (q.num_quads,), f"{name}: a scaled Jacobian per quad")
        before = q.quality["min_scaled_jacobian"]
        r = q.smooth(20)
        check(r["ran"] and r["inverted_after"] == 0,
              f"{name}: smooth(20) {before:.3f} -> {r['min_scaled_jacobian_after']:.3f}, no folds")
        check(abs(q.quality["min_scaled_jacobian"] - r["min_scaled_jacobian_after"]) < 1e-12,
              f"{name}: the smoothed mesh is the one the object holds")
    print("all checks passed")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/meshes/singlemat/geom001.obj")
