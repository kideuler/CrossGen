# Tier-1 experiments, E1 to E5

The experiments of Sec. 6 of the IMR 2027 submission, as five programs. Each one
prints a table you can read while it runs and writes CSVs a plot script or the
supplement can be built from; each exits non-zero if a check it makes fails, so
they double as regression tests. E6 (the Laghos bubble run) lives outside this
repository and is not here.

## 2026-09-03 (latest): the mechanism set, `data/meshes/mechanism/`

Eleven domains built to test *where the two discretisations must differ*, built
by `MakeDomains` from `paper::mechanismDomains` and written as `.obj` so every
experiment can run on them with `--set mechanism` (E5 now takes `--set` too).

**They are deliberately not in `singlemat/` or `multimat/`.** Those are the
corpus the paper's aggregate tables average over, and a corpus chosen after
seeing the results is not evidence. What makes this set legitimate rather than
cherry-picked is that the prediction is derived from the code before the run,
each family sweeps one parameter, and each family contains the setting at which
the prediction says the two methods *agree*.

**The prediction.** `CrossField` prescribes its Dirichlet data per *vertex*, by
rounding a boundary corner's interior angle into one of four quarter-turn
classes (`CrossField.cxx:114-127`). That rounding is exact at 90, 180 and 270
degrees and wrong by up to 22.5 degrees half way between. `DualMBO` pins per
*edge*, to that edge's own tangent, and rounds nothing. Remark 4.1 is the same
argument at a material junction. So the two should separate wherever a corner or
a junction sector is far from a multiple of a right angle, and agree where it is
not.

| family | sweeps | control |
|---|---|---|
| `comb` | nothing: all angles 90 or 270 | this is the control |
| `star5t{90,70,55,45}` | tip angle, corner deviation 18 to 45 deg | none; `comb` is |
| `laminate{90,75,60,45}` | junction obliquity 0 to 45 deg, *square boundary* | `laminate90` |
| `grain{7,19}` | Voronoi polycrystal, junctions generically oblique | none |

`laminate` is the sharp one: the domain is the unit square, so every boundary
corner is a right angle and the only oblique thing in it is where an interface
lands on the boundary. It isolates the junction mechanism from the corner one.

**What was measured** (`results/mech/`, two-point weight):

* B1's 95th-percentile interface misalignment **equals the junction obliquity**,
  to 2e-3 degrees, on all four laminates (0, 15, 30, 45). Ours is 1e-14 on all
  four. On the polycrystals B1 is 39 degrees and ours 1.3e-14.
* On boundary corners: at 18 and 38 degrees of corner deviation ours is 1e-14
  and B1 is 4.9 and 9.7 degrees. At 45 degrees **ours is worse** (40.5 against
  22.5) -- the corner-face degeneracy of Sec. 4.4, which the set was built to
  catch and does.
* Controls tie exactly: `comb` and `laminate90` give both methods the same
  layout, 272 and 400 quads at minimum scaled Jacobian 1.000.
* Layout: mostly a wash, as on the corpus, with two exceptions in our favour --
  `star5t55` (0.714 against 0.508) and `laminate60` (1 inverted quad against
  52) -- and one on the polycrystal, `grain7`, where neither layout validates
  but ours has 14 inverted quads against 147.
* `grain7`'s layout does not complete for either field, which is not a surprise:
  E5(b) put the quantization stages' threshold at about 110 degrees of sector
  and Voronoi triple junctions sit at about 120.

## 2026-09-03 (later): five additions, and one defect in the shipped pipeline

Everything below the next heading still holds. What changed after it:

* **The pipeline was not running the method the paper describes.** `TORSION`'s
  Stage 0 constructed `DualMBO` with the default `MinHeight` weight and took a
  single `tau = D^2/10` -- neither the two-point weight of Sec. 4.3 nor the
  continuation of Sec. 4.5. E4/E5 were unaffected (they substitute an external
  field and Stage 0 is skipped), but anything produced by running the pipeline
  itself was, and that includes the hydro mesh. `TORSION::Options` now carries
  `dualMBOWeight` (default `Orthogonal`) and `dualMBOTauContinuation` /
  `dualMBOTauRatio` / `dualMBOTauFloorEdges` / `dualMBOTauLevelSteps`, and
  `TestTORSION` takes `--weight min|harm|orth`, `--no-continuation`,
  `--tau-ratio`, `--tau-floor`. On `bubbles` with disk templates the fix takes
  the worst scaled Jacobian from 0.748 to 0.818.

* **The SparseLU fallback is gone.** Nine of the 24 single-material models (and
  `bubbles`, `rocket`) carry vertices no triangle references -- `geom031` has
  309 -- which leave an empty row in the P1 mass and stiffness matrices, so
  `A` is structurally singular and Eigen refuses the factorisation.
  `CrossField::initialize` now gives such a vertex a unit mass entry and pins
  it like a boundary vertex. The field is unchanged (one B1 run differs by one
  iteration over the whole corpus, from BiCGSTAB's tolerance no longer being in
  the loop); every run in `results/` now factorises directly.

* **E1(f), the graded-mesh sweep** (`E1_graded.csv`): `gradedDiskDomain(radius,
  rings, growth)` builds the disk as `rings` concentric rings whose radial
  spacing grows geometrically, as the constrained Delaunay triangulation of
  exactly those points (`TriangleMesher2D::Options::just_delaunay`). Both
  weights run on the *same* mesh at each grading, and the sweep reports the
  affine-field residual beside the cone-placement error against the closed-form
  answer. This is the experiment the paper needed: on the near-uniform corpus
  the two weights give nearly the same field, and the difference only becomes
  visible once the mesh is graded (min-height 21.3 degrees of cone error at
  g = 1.2 against the two-point weight's 2.4).

* **E5(b), the oblique-junction sweep** (`E5_oblique_sweep.csv`, via
  `--oblique-sweep 90,100,105,110,115,120`): `obliqueJunctionDomain` now takes a
  sector angle. The field is exact at every angle; the layout completes a valid
  mesh up to 110 degrees and stops short at 115 and 120. That is what the
  `oblique` failure is -- a threshold in Stages 6-8, not a defect of the field.

* **Fair-comparison and provenance switches.** `E3_Benchmark` takes
  `--b1-continuation` (give the baseline our tau ladder; results in
  `results/orth_b1cont/`) and `--b1-seeds N` (re-run B1 from N random starts and
  report how much of its singularity set is the seed -- on this corpus, none of
  it). `E4_Layout` and `E5_MultiMaterial` take `--layout-opt NAME=VALUE` to
  sweep a `TORSION::Options` field by name, and print the pipeline's own stage
  messages for any layout that did not validate. `LayoutResult` now carries
  `externalFieldUsed`: with `--disk-templates` on a domain with a circular
  inclusion, Stage 0c changes the triangle count, the supplied field is refused
  and the pipeline solves its own -- such a row is *not* a comparison of the two
  methods and E5 now says so. `LayoutResult` also carries Stage 4's integration
  fit residual, which is the one field-quality number the layout stage measures
  for us.

Results directories after this work:

```
paper_tests/results/orth/           two-point weight -- the method as shipped
paper_tests/results/min/            min-height weight -- the ablation rows
paper_tests/results/orth_b1cont/    two-point, with B1 given the tau ladder
paper_tests/results/orth_templates/ two-point, E5 with the disk templates on
```

## 2026-09-03: the penalty weight is a switch, and the results are per weight

Every experiment now takes `--weight min|harm|orth` (the DualMBO penalty weight
of `DualMBO::PenaltyWeight`; `min` is the interior-penalty default, `orth` the
two-point / circumcentric-dual weight the paper ships). It is a process-wide
default (`paper::setDefaultPenaltyWeight`) so that every `MethodOptions` an
experiment constructs carries it, and each log's header line says which weight
ran. The results the paper is written from live in

```
paper_tests/results/orth/   the two-point weight  -- the method as shipped
paper_tests/results/min/    the min-height weight -- the ablation rows
```

each with its own `E*.log`, CSVs, `vtk/` and `obj/`. `paper_tests/results/*.csv`
at the top level are an older min-height run and are superseded.

Two additions to the harness the paper's Sec. 6 depends on:

* **E1(e), consistency** (`E1_consistency.csv`): an affine field sampled at the
  circumcenters, `rho = |M^-1 K phi| h_mean / gamma` read on every triangle whose
  three edges are interior, for all three weights, over both corpora; plus the
  per-mesh counts of non-Delaunay and floor-clamped edges. The two-point residual
  is asserted away from floor-clamped triangles (Prop. 4 is a statement about the
  unclamped operator) and reported on them. Uses the new read-only accessors
  `DualMBO::stiffnessMatrix()` / `massMatrix()`.
* **E3's constrained energy** (`E3_<set>_energy.csv`, and the `energy_twopoint`,
  `energy_boundary` columns of the per-model CSV): the interior energy in the
  two-point functional (`metrics::edgeWeightsTwoPoint`, neutral between
  operators) and the one-sided constraint term (`metrics::boundaryEnergy`), so
  that constraint violation is charged. This is the table that reverses the
  interior-only energy comparison.

Numbers quoted further down this file are from the earlier top-level min-height
run and are kept for the record; the paper's numbers are the ones in
`results/orth/` and `results/min/`.

```
paper_tests/
  common/       the shared machinery, so that every experiment measures the
                same quantity the same way
    Domains     the canonical domains of E1/E2 and the authored junction
                domains of E5, built in memory
    Methods     the three methods behind one interface, and the documented
                B1 conversion
    Metrics     Sec. 5's definitions, once
    Layout      E4/E5's end of the pipeline: a field in, a quad mesh out,
                through TORSION
    Corpus      data/meshes enumeration
    Report      CSV writer and console tables
  E1_Verification.cxx
  E2_Canonical.cxx
  E3_Benchmark.cxx
  E4_Layout.cxx
  E5_MultiMaterial.cxx
```

## Running them

```
cmake --build build --target PaperE1_Verification PaperE2_Canonical \
      PaperE3_Benchmark PaperE4_Layout PaperE5_MultiMaterial
cd build/paper_tests

./E1_Verification --out results --levels 5      # ~1 min at 4 levels, ~10 at 5
./E2_Canonical    --out results --vtk results/vtk
./E3_Benchmark    --out results                 # ~30 s, 24 models x 3 methods
./E4_Layout       --out results --obj results/obj    # the long one: the pipeline
./E5_MultiMaterial --out results --obj results/obj --vtk results/vtk
```

or `make paper` for all five. `--help` on any of them lists its switches. E1 and
E2 are also ctest cases (`ctest -R Paper_`); E3 to E5 are not, because a full
sweep is minutes to hours and belongs in a run of its own.

## The three methods

Every experiment compares the same three, and all three produce the same object:
one unit spin-4 value `u = e^{4i theta}` per triangle.

| | what it is | DOFs | singularities |
|---|---|---|---|
| **DualMBO** | ours: p=0 dual-mesh MBO (`src/dualmbo`) | faces | vertices, natively |
| **B1** | P1-MBO (`src/crossfield/CrossField`) + the documented conversion | vertices | faces, then resampled |
| **B2** | polyvector field (`src/polyvector/PolyVectors`) | faces | vertices |

**The B1 conversion** is `paper::convertP1ToFaces` and is one function so that the
paper can describe it in half a column: the face value is the spin-4 average of
the three corner values, renormalised; per-edge matchings are the quarter turn
between two face values; a singular face is mapped to the vertex of that face
nearest its centroid. What it costs is measured, not asserted --
`metrics::ConversionReport`.

**Fair-comparison protocol.** Energies are always the face functional of Sec. 5
at a fixed evaluation weight `gamma = 10`, whatever gamma a solver ran at, and
interface edges are excluded from it because the DualMBO operator carries no
coupling there. B1 uses a seeded random start and E2 reports mean and sigma over
ten seeds. Timings are single-threaded Eigen SparseLU on one machine.

## Findings the experiments produced that the paper should use

These came out of running the code and are not in the outline.

1. **gamma and tau are one parameter, not two.** For p=0 the stiffness matrix and
   the boundary right-hand side are both exactly proportional to gamma, and both
   enter the MBO step only through `A = M + tau K` and `r = M u + tau b`. So the
   iteration depends on the two only through the product `gamma * tau`.
   E1(a/b) verifies it to 1e-13 on four domains. This belongs in Sec. 4.5 and it
   collapses two sweeps into one.

2. **The shipped stopping criterion stops before the energy has settled.**
   `error < 2 N 1e-5` fires after two or three steps at `tau = D^2/10`, at an
   energy about 5% above the fixed point's on the unit disk. The singularity
   *structure* is already right there, which is why the pipeline is happy with
   it, but an energy comparison is not. The field-quality experiments (E2, E3, E5
   part 1) therefore run at 1e-9 and say so; the pipeline experiments (E4, E5
   part 3) leave it at the shipped value, because what they measure is the
   pipeline.

3. **Prop. 2 has to be stated in the heat content, not in a Dirichlet energy.**
   MBO is monotone in the Esedoglu-Otto functional; the Dirichlet-type energy is
   not what threshold dynamics decreases. Run to convergence the two agree on
   this corpus (E1(d) measures the worst rise at 4e-16 relative), but the
   statement should be made about the functional the argument is about.

4. **Poincare-Hopf needs the field's own boundary index, and the corpus needed a
   fix to compute it.** Two traps, both now handled in `metrics::singularities`:
   several models in `data/meshes` carry vertices no face references (geom031 has
   309), and counting them puts chi out by exactly that many; and the boundary
   fan cannot be read off `VertexTriangleCSR`, which sorts a star by centroid
   angle about a branch cut unrelated to where the boundary is -- it has to be
   walked by edge adjacency. With both fixed, **Prop. 3 holds on all 24
   single-material models for DualMBO and B2, and for B1 in the field form on all
   24 and in the geometric form on 20 of 24** (E3_singlemat_aggregate.csv).

5. **B1 has no interface-alignment mechanism at all.** `CrossField`'s only
   Dirichlet data is dS. Its interface alignment error on a multi-material domain
   is therefore the absence of the constraint rather than a worse answer to the
   same question, and its lower energy on those domains is the energy of an
   unconstrained problem. E5 says this in its output; the paper must say it too
   or a reviewer will call the comparison unfair.

6. **A junction is only a conflict when it is oblique, and E5's headline numbers
   fall straight out of that.** A cross is invariant
   under a quarter turn, so a junction all of whose sectors are multiples of 90
   degrees is satisfied by one cross and a vertex-based field has no difficulty
   there. E5 reports a per-junction **obliquity** so that the proposition of
   Sec. 4.4 is read on the rows where it applies. Measured:

   | domain | junction obliquity | DualMBO junction residual | B1 junction residual |
   |---|---|---|---|
   | `tjunction` (control) | 0 deg | 7e-15 | 1.2e-05 |
   | `quadruple` (control) | 0 deg | 7e-15 | 2.7e-07 |
   | `oblique` | 30 deg | 1.2e-14 | **35.2 deg** |
   | `geom004` (corpus) | 30 deg | 15 deg* | **38.0 deg** |

   Over all 14 domains and 61 junctions, DualMBO's residual maxes at 16.9 deg
   (mean 1.58 deg, 12 junctions over 1 deg) against B1's 38.0 deg (mean 6.06 deg,
   33 over 1 deg) (E5.log, "Part 2 summary").

   \* DualMBO's 15 degrees on geom004 is two edge-sides out of 210, at triangles
   carrying two non-parallel interface edges; its 95th percentile is 7e-15. That
   is the same corner-averaging limitation Sec. 7 lists for dS, and E5 counts it
   rather than hiding it behind a percentile.

   And R1, the number the application audience cares about: **the mixed-zone
   count is 0 on every one of the 14 domains, for both fields** -- including the
   `mesh(partial)` runs (`bubbles`, `inclusion`, and others) that used to be the
   exceptions. The one domain R1 is not measurable on is `oblique`/DualMBO,
   where the layout produced no quads at all. This is a real improvement, not
   just a smaller corpus: on the previous corpus `bubbles`, `icf`, and
   `inclusion` all had nonzero mixed-zone counts; `icf` is gone from the trimmed
   corpus, but `bubbles` and `inclusion` are still in it and now measure 0,
   consistent with Stage 2's cut now being routed around the material interface
   network rather than through it.

7. **E4 does not support "equal-or-better quad layouts", and the paper must not
   claim it.** On the 24 single-material models, both fields reach a mesh on all
   24; layouts are valid on 23 for each field -- this used to favour DualMBO
   (34 against 33 of 35) but on the trimmed corpus it's a tie. The paired
   scoreboard on final element quality still goes mostly the other way:

   | metric | DualMBO better | tie | B1 better |
   |---|---|---|---|
   | min scaled Jacobian | 4 | 17 | 3 |
   | mean scaled Jacobian | 3 | 15 | 6 |
   | inverted quads | 1 | 22 | 1 |
   | irregular interior | 1 | 19 | 4 |
   | min zone dimension | 5 | 11 | 8 |
   | 1% zone dimension | 6 | 13 | 5 |
   | max angle deviation | 6 | 13 | 5 |

   DualMBO still carries more singularities than B1 on average (6.88 against
   5.92, `E3_singlemat_aggregate.csv`), but on this corpus B1's layouts now have
   *more* total patches (1058 against DualMBO's 893 over the same 24 models,
   `E4.log`) even though DualMBO has more interior cones in total (165 against
   142). **This reverses the old "more cones -> more patches" mechanism claim
   from the pre-trim corpus and needs a look before it goes in the paper** --
   the paired scoreboard and mesh/valid counts above are read straight off the
   fresh run, but no explanation for the patches reversal has been verified
   here. B1's geometric Poincare-Hopf check still fails on some models (4 of
   24, same absolute count as before, now a larger fraction: 20/24 pass).

   So the claim E4 supports is about **robustness and topological consistency**,
   not element quality:

   * every model reaches a mesh on both fields (24/24); valid layouts tie at
     23/24 each (no longer a DualMBO edge on this corpus);
   * the field's topological content is what the quantizer receives, with no
     conversion in between -- against 16 lost edge matchings over 9 models and
     8 singularities created over 3 models for B1 (E3);
   * Prop. 3 holds for DualMBO on every model in both forms, and for B1 in the
     geometric form on only 20 of 24.

   Sec. 1's contribution bullet 5 should be rewritten to that, and to drop the
   "DualMBO reaches one more valid layout" framing now that the corpus ties. If
   element quality is wanted as a claim, the thing to investigate first is
   DualMBO's extra cones -- E2 shows the corner-averaging Dirichlet data is
   degenerate at a 45-degree corner (45 degrees of alignment error there, up
   from a previously-measured 41.8), which is a plausible source of cones the
   geometry does not need.

## What each experiment writes

| file | contents |
|---|---|
| `E1_gamma.csv`, `E1_tau.csv` | the two sweeps |
| `E1_product.csv` | the gamma*tau identity check |
| `E1_refine.csv` | h-refinement, energy and cone positions |
| `E1_convergence.csv` | per-iteration increment and both energies |
| `E1_consistency.csv` | E1(e): affine-field residual per mesh and weight, degenerate-edge counts |
| `E2_canonical.csv` | the canonical-domain table, per method |
| `E2_robustness.csv` | B1's seeds and the jittered remeshes |
| `E3_<set>_permodel.csv` | every model, every method |
| `E3_<set>_aggregate.csv` | the aggregate row per method |
| `E3_<set>_conversion.csv` | the conversion-cost totals |
| `E3_<set>_energy.csv` | interior (min-height and two-point) and boundary energies summed per method |
| `E4_<set>_permodel.csv` | the pipeline outcome per model and field |
| `E4_<set>_aggregate.csv` | success rates and quality aggregates |
| `E5_field.csv` | interface alignment and junction residuals |
| `E5_layout.csv` | mixed zones, irregular vertices per material, quality |
| `E5a_bc_ablation.csv` | hard-pinned against weak Nitsche |

`--vtk DIR` writes a legacy VTK per field (cross axes as cell vectors, index as
point data) and `--obj DIR` writes the final quad meshes, for the figures.

## The two hooks this needed in `src/`

Both are small and both are documented where they live.

* `DualMBO::setTauScale` / `getTau` / `getGamma` (E1's sweeps) and
  `DualMBO::setHardBoundaryConditions` (E5(a)'s ablation -- the weak Nitsche form is
  the same assembly with the Dirichlet row elimination left out).
* `TORSION::Options::externalField`: a per-face field to use instead of running
  Stage 0. Every stage downstream reads the field through `DualMBO::u_k_prev` and
  nothing else, so substituting the vector is the whole of substituting the
  field, and that is what makes E4 a controlled experiment rather than a
  comparison of two pipelines.
