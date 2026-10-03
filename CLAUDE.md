# CrossGen

2D cross fields, and quad layouts / block decompositions built from them, for
planar multi-material domains. C++17, Eigen (submodule), OpenCASCADE, Gmsh,
optional Qt6 viewer. It backs an IMR 2027 paper: a p=0 dual-mesh MBO cross field
(`src/dualmbo`) feeding the Shepherd–Gu–Hughes 2022 quad-layout pipeline.

## Rules

- **Read-only git only.** `status`/`diff`/`log`/`show` are fine. Never `add`, `commit`,
  `push`, `stash`, `checkout`, `reset`, `restore` or `clean`. The user manages history.
  If you need a "before" baseline, add a runtime flag or env-var override instead of stashing.
- **OpenCASCADE is confined to `src/geom`.** No OCC header, and no `geom/detail/`
  include, anywhere else. The build enforces it (OCC is linked PRIVATE to
  `CrossGenGeom`) and so does ctest `Geom_OpenCascadeConfinedToGeom`. New geometry
  goes into `src/geom`, behind public headers that name no OCC type.
- **Ignore `paper_tests/`.** Don't read it, run it, log findings in its README, or wire
  new options into it. Its `Paper_` ctests are not a signal either. Jacob looks after it.
- Put throwaway drivers, corpus sweeps and dumps in the session scratchpad, not in the repo.
- Match the house style: `.hxx`/`.cxx`, `#ifndef __NAME_HXX__` guards, include paths
  relative to `src/` (`"MERIDIAN/Immersion.hxx"`). Comments are long-form prose that
  explains *why* and cites paper sections (e.g. "Sec. 3.3", "Q1–Q5", "Definition 2.1").

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release          # existing build/ is Release, viewer ON
cmake --build build -j8                                 # or: --target TestMERIDIAN TestTORSION -j8
```

- `build.sh` runs `rm -rf build` and does a full rebuild. `make_meshes.sh` also calls
  `clean.sh`, which deletes the `.msh`/`.obj` meshes. Don't run either unless asked.
- Options: `BUILD_OPENGL_VIEWER` (Qt6 `Viewer`), `CROSSGEN_ENABLE_OPENMP` (TMOP; finds Homebrew
  libomp), `BUILD_MESH2GMSH`, `BUILD_PYTHON_MODULE` (ON: the `crossgen` extension module in
  `build/python`, built for the `python3` on PATH; use it with `PYTHONPATH=build/python`, or configure once with `-DCROSSGEN_PYTHON_DEV_PTH=ON`, which writes `crossgen-dev.pth` into that Python's site-packages so no `PYTHONPATH` is needed).
- Libraries: `CrossGenGeom` (src/geom, OCC) and `PolyVector` (all the rest, links
  CrossGenGeom PUBLIC). Every executable is a thin driver in `src/utils/`.
- `data/meshes/**` is copied into `build/data/meshes` **at configure time**. Re-run
  `cmake` after adding or regenerating a mesh.

## Test

- `cd build && ctest` (or `ctest -R MERIDIAN`, `-R TORSION`, `-R ZIPLINE`, `-R Geom`, `-R ShapeDNA`, `-R Python`;
  `ctest -E Paper_` leaves out the paper_tests cases, which are ignored).
  A project PreToolUse hook (`.claude/hooks/filter-test-output.sh`) cuts `ctest`/`make test`
  output down to the FAIL/ERROR lines. The end-to-end cases have 2400 s timeouts; the
  bubbles template cases are the slow ones.
- Self-tests: `TestMERIDIAN --selftest`, `TestTORSION --selftest`, `TestTMOP --selftest`, `TestGeom`, `TestZIPLINE --selftest`, `TestShapeDNA --selftest`.
- A single model: `./build/TestMERIDIAN data/meshes/singlemat/geom012.obj [flags]`
  (TestTORSION works the same way). The flag parser is in `src/utils/TestMERIDIAN.cxx`
  around line 330. Useful flags: `--disk-templates --near-miss 0.04` (bubbles),
  `--tmop N`, `--cones N`, `--curves N`, `--step/--brep` (write the output),
  `--sep`/`--sep-uv`/`--psi`/`--layout` (OBJ dumps). TestTORSION adds
  `--weight min|harm|orth`, `--no-continuation`, `--ref-test`.
- Whole corpus: `make meridian` / `make torsion` / `make dualmbo` run from `build/`
  and don't fail on a bad model. For comparisons, write a scratch driver and run it with `xargs -P`.
- `TopoDiag <mesh>` is the user's own diagnostic for a single model.

## Layout

| Path | What |
|---|---|
| `src/mesh` | `Mesh` (triangle mesh, material tags), `QuadMesh` + MFEM `.mesh` writer, `TMOP` smoother, `Pillow` (a layer of quads along dS/an interface wherever an element spans a flat feature corner, which no smoothing can repair; the pipelines' Stage 12 runs it before TMOP), `BoundaryFeatures` (corners by the 20° rule, holes, T3's `minimumDefect` that `BlockDecomposition::regularity()` grades against), `FeatureFrame` (the harmonic extension of the walls' and interfaces' tangent cross, which `alignmentQuality()` grades blocks and meshes against) |
| `src/dualmbo` | **DualMBO**, our method: p=0 dual-mesh MBO. Crosses on faces, singularities on vertices |
| `src/crossfield`, `src/polyvector` | Baselines B1 (P1-MBO, per vertex) and B2 (polyvectors) |
| `src/MERIDIAN` | **Pipeline A**, Shepherd 2022 Stages 0b–11: Interfaces → ConeSingularities → ConeCut → RicciFlow → Immersion → SubdomainLabels → LayoutEnergy → Separatrices → Arrangement → SplineFit → QuadMesh, plus DiskTemplate (0c/11). The stage map is in the header comment of `MERIDIAN.hxx` |
| `src/TORSION` | **Pipeline B**: the same stages, but Stages 3–4 are swapped for integrating the DualMBO field (ConeMetric, FieldFrames, FieldIntegration, TutteEmbedding). Both pipelines meet at `Immersion`, and every stage after that is shared. On a multi-material model it lays each material region out on its own by default (`MaterialLayout`, `Options::perMaterial`; `--whole-model` / `per_material=False` for the old route): Stages 1–8 per region, matched across the interfaces, glued into one `Arrangement` for Stages 9–11, with the whole-model run as fallback. The viewer's TORSION mode does the same ('w' switches before the cone phase), on TORSION's own tau-continued field (`TORSION::prepareFieldLevel`/`fieldTauLadder`; MERIDIAN mode keeps MERIDIAN's single-tau one) |
| `src/geom` | Splines, polylines, Coons, B-rep (`Vertex/Edge/Face/Shape`), STEP/BREP output, all on top of OCC |
| `src/ZIPLINE` | **ZIPLINE**, Viertel IMR19 (class `ZIPLINE` in `ZIPLINE.hxx` runs every stage; viewer mode 2 and driver `TestZIPLINE` both go through it): separatrices → T-layout → Sec. 4 simplification + Sec. 12 stem extension → `LayoutBlocks` (the shared BlockDecomposition on src/geom splines; faces with a T-junction are not blocks, so read it by coverage) → BlockQuadMesh/TMOP. Honours material interfaces through MERIDIAN's Stage 0b `Interfaces` and an interface-aligned `CrossField`. `docs/viertel_2019.md` Sec. 14 maps the spec to the code |
| `src/quantization`, `src/medialaxis`, `src/UMBER`, `src/OASIS`, `src/Parameterization` | Older approaches (QGP, medial axis, polysquare, spectral, seam cuts). Rarely touched now |
| `src/viewer` | Qt6 `Viewer`. `run_meshes.sh` opens every mesh in turn; it has figure/SVG export |
| `src/python` | **`crossgen`**, the CPython extension module (plain `Python.h`, no pybind): `crossgen.load(obj)` → `Mesh`; `m.zipline()/umber()/meridian()/torsion()/atlas()` → `BlockDecomposition`; `b.mesh(h)` → `QuadMesh`; `q.smooth(n)` (TMOP); `m.shape_dna()` → ndarray; `m.boundary_features()` → dict; `b.valid`; `b.regularity()/angle_quality()/chord_quality()/alignment_quality()` (the C++ `BlockDecomposition` methods, qualities in [0, 1], higher better, NaN when `b.valid` is False; `docs/block_decomposition_metrics.md` Secs. 6-7; `alignment_quality` is how the blocks' grids follow the walls, against `m.feature_frame()`, and `q.alignment_quality()` reads it on a mesh). `py/build_dataset.py` writes the selector's training CSVs from these (the `shape_dna(normalization='weyl_ratio')` spectrum only with `--eigs K`: it made the selector worse and was a third of the build), each method averaged over `--copies` exact symmetries of the mesh, plus ZIPLINE's first run on its own as `probe_zipline_*`; a rerun reuses `cache.jsonl` and runs only missing tasks. `py/train_classifier.py` trains on them (a validity head and a quality head for one metric, `--metric`, default `alignment_quality`) for the probe rule: run ZIPLINE, keep it if it ranks first, else run the top other method and keep the better. Keywords are the C++ `Options` fields in snake_case, from one table per struct (`PyOptions.hxx`); `crossgen.options(method, stage)` lists them. Each method meshes through its own mesher (`Methods.hxx` says why). UMBER has no driver class, so `MethodUMBER.cxx` chains its stages the way `TestUMBER` does, and `MethodLayout.cxx` rebuilds Stages 10–11 the way `MERIDIAN::run()` does. Change those two files if the driver or `run()` changes. ctest `Python_Module` |
| `paper_tests/` | E1–E5 paper experiments. **Ignored** (see Rules): not read, run, logged in or extended |
| `data/geometry/{singlemat,multimat}` | `.geo` sources. `Mesh2Dgmsh` turns them into `data/meshes/<kind>/*.obj` (gitignored) |
| `data/meshes/{singlemat,multimat,mechanism}` | Corpus: 24 + 11 + mechanism set (built by `MakeDomains`, kept separate from the corpus on purpose) |
| `data/geometry_extra/` | Models dropped from the corpus. Not used |
| `src/ShapeDNA` | **Shape-DNA** (Reuter et al. 2006, `docs/shape_dna.md`): the first n Laplacian eigenvalues of the whole domain (materials ignored), normalised for scale. P1–P3 Lagrange FEM with exact reference integrals → our own nested-dissection multifrontal Cholesky → our own shift-invert block Lanczos (thick restart, B-inner product). All OpenMP, and bit-identical for any thread count. Driver `TestShapeDNA <mesh.obj>...` prints each DNA and the pairwise distances. P3/50 eigenvalues on the largest corpus model takes ~1.8 s |
| `docs/` | Design notes: `ricci_flow_pipeline.md` (A), `cf_flow_pipeline.md` (B), `multimaterial.md`, `tmop_node_mobility.md`, `shape_dna.md`, `shepherd2022.pdf` |

## Before changing a pipeline stage

The memory index has a note for each area (stages 6–10, ZIPLINE, multimat, disk
templates, TMOP, SIPG/DualMBO, geom/OCC). Read the note for the area before
changing it. Most of them record a trap that already cost a debugging session.
