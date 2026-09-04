// E1 -- Verification and parameter study (Sec. 6, 0.6 page).
//
//   (a) gamma sweep      {1, 5, 10, 20, 50, 100}: final energy, #singularities,
//                        MBO iterations. The claim under test is a plateau of
//                        insensitivity wide enough that gamma = 10 is a choice
//                        and not a tuning.
//   (b) tau sweep        x{0.1, 0.5, 1, 2, 10} about the D^2/10 heuristic:
//                        singularity count against over- and under-diffusion.
//   (ab) the product     gamma and tau are **not two parameters**. For p = 0 the
//                        stiffness matrix and the boundary right-hand side are
//                        both exactly proportional to gamma -- K = gamma K_0 and
//                        b = gamma b_0 -- and both enter the step only through
//                        A = M + tau K and r = M u + tau b, so the whole
//                        iteration depends on gamma and tau through the single
//                        product gamma*tau. This part verifies that to machine
//                        precision, which turns (a) and (b) from two sweeps into
//                        one sweep read twice and is worth a sentence in Sec. 4.5.
//   (c) h-refinement     the common energy and the cone positions against h on
//                        the disk, whose answer is known -- four +1/4 cones at
//                        45 + k*90 degrees from the centre, at a radius that
//                        goes to zero with h.
//   (d) MBO convergence  the increment ||u^{k+1} - u^k|| and two energies per
//                        iteration, over a fixed budget rather than to the
//                        convergence test. This is where Prop. 2 is checked,
//                        and it is reported as a measurement: which functional
//                        actually decreased, and by how much any rise violated
//                        it. A Dirichlet-type energy is *not* the functional
//                        MBO is monotone in -- threshold dynamics decreases the
//                        heat content of Esedoglu and Otto, not the Dirichlet
//                        energy -- so a small rise here is a fact about how
//                        Prop. 2 has to be stated, not a bug.
//   (e) consistency      at p = 0 the edge weight kappa_e *is* the discrete
//                        Laplacian, and a second-order operator has to
//                        annihilate affine fields. An affine phi is sampled at
//                        the circumcenters of every mesh in the corpus and
//                        M^-1 K phi is read on the triangles whose three edges
//                        are interior, for each of the three weights DualMBO
//                        can assemble. The two-point (circumcentric dual)
//                        weight gives rounding by Prop. 4; the interior-penalty
//                        weight does not, and this part measures by how much.
//                        It also counts the degenerate edges -- cocircular
//                        pairs, where the dual edge has no length -- that the
//                        floor kMinDualDistance exists for.
//
// The comparison energies are measured with the *same* functional -- the
// face-based one of Sec. 5 at the fixed evaluation weight gamma_eval = 10 --
// whatever gamma the solver ran at. Comparing a solve at gamma = 100 against
// its own energy would compare two different functionals and would say nothing.
//
// `--weight min|harm|orth` selects the penalty weight parts (a) to (d) run
// with; (e) always measures all three.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "common/Corpus.hxx"
#include "common/Domains.hxx"
#include "common/Methods.hxx"
#include "common/Metrics.hxx"
#include "common/Report.hxx"

using namespace paper;

namespace {

// The weight every energy in this file is evaluated at, whatever the solver ran
// at. Sec. 5: "DualMBO weights fixed for all methods".
constexpr double kGammaEval = 10.0;

struct Solved {
    FieldRun run;
    metrics::SingularityReport sing;
    double energy = 0.0;
};

Solved solve(const std::shared_ptr<Mesh> &m, const MethodOptions &o) {
    Solved s;
    // E1 is the study of the *single-tau* operator, and every part of it means
    // that: (a) and (b) sweep one tau, the product identity of (a/b) is a
    // statement about one step -- A = M + (gamma tau) K_0 -- and (c) refines at
    // fixed tau. The tau-continuation the other experiments run would make all
    // three measure something else, and the product identity would simply be
    // false under it, since the continuation's floor depends on gamma and on
    // the mesh separately rather than on gamma*tau. So it is off here, and what
    // the continuation buys is measured on its own in E1(e).
    MethodOptions single = o;
    single.tauContinuation = false;
    s.run = runDualMBO(m, single);
    if (!s.run.ok) return s;
    s.sing = metrics::singularities(*m, s.run.u);
    s.energy = metrics::commonEnergy(*m, metrics::edgeWeights(*m, kGammaEval), s.run.u,
                                     interfaceEdgeFlags(*m));
    return s;
}

// The mesh sizes E1(c) refines through. Chosen so that the coarsest is about a
// thousand triangles and each level is roughly four times the last.
std::vector<double> refinementSizes(int levels) {
    std::vector<double> h;
    double x = 0.06;
    for (int i = 0; i < levels; ++i) { h.push_back(x); x *= 0.5; }
    return h;
}

// The circumcenter, exactly as DualMBO::initialize computes it.
Point circumcenterOf(const Point &pa, const Point &pb, const Point &pc) {
    const double bx = pb[0] - pa[0], by = pb[1] - pa[1];
    const double cx = pc[0] - pa[0], cy = pc[1] - pa[1];
    const double den = 2.0 * (bx * cy - by * cx);
    const double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    Point cc = pa;
    if (std::fabs(den) > 1e-300) {
        cc[0] = pa[0] + (cy * b2 - by * c2) / den;
        cc[1] = pa[1] + (bx * c2 - cx * b2) / den;
    }
    return cc;
}

// Everything about a mesh's circumcentric dual that the consistency test and
// the grading sweep both need: where the cell centres are, how long each dual
// edge is against the element size, which dual edges the floor clamps, how
// unequal the two elements on an edge are, and which triangles have three
// interior edges and so can be read without a boundary row in the way.
struct DualStats {
    std::vector<Point> circ, cent;
    std::vector<double> area;
    std::vector<char> clampedEdge;   // per edge
    std::vector<int> tested;         // triangles with three interior edges
    std::vector<char> clamped;       // parallel to `tested`
    std::vector<double> dualRatios;  // d_e / d_ref, per interior edge
    int interiorEdges = 0, nonDelaunay = 0, belowFloor = 0, clampedTriangles = 0;
    double maxAreaRatio = 1.0, medianDualRatio = 0.0, hMean = 1.0;
};

DualStats dualStats(const Mesh &m) {
    DualStats d;
    const int NT = static_cast<int>(m.triangles.size());
    d.circ.resize(NT);
    d.cent.resize(NT);
    d.area.resize(NT);
    for (int t = 0; t < NT; ++t) {
        const Triangle &tri = m.triangles[t];
        const Point &pa = m.vertices[tri[0]];
        const Point &pb = m.vertices[tri[1]];
        const Point &pc = m.vertices[tri[2]];
        d.circ[t] = circumcenterOf(pa, pb, pc);
        d.cent[t][0] = (pa[0] + pb[0] + pc[0]) / 3.0;
        d.cent[t][1] = (pa[1] + pb[1] + pc[1]) / 3.0;
        d.area[t] = 0.5 * std::fabs((pb[0] - pa[0]) * (pc[1] - pa[1])
                                  - (pc[0] - pa[0]) * (pb[1] - pa[1]));
    }

    d.clampedEdge.assign(m.edges.size(), 0);
    double hSum = 0.0;
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        const Point &ea = m.vertices[m.edges[e][0]];
        const Point &eb = m.vertices[m.edges[e][1]];
        const double len = std::hypot(eb[0] - ea[0], eb[1] - ea[1]);
        hSum += len;
        if (m.isBoundaryEdge[e] || len < 1e-14) continue;
        const int ti = m.edgeTriangles[e][0], tj = m.edgeTriangles[e][1];
        if (ti < 0 || tj < 0) continue;
        ++d.interiorEdges;
        const double hi = 2.0 * d.area[ti] / len, hj = 2.0 * d.area[tj] / len;
        double nx = -(eb[1] - ea[1]) / len, ny = (eb[0] - ea[0]) / len;
        if ((d.cent[tj][0] - d.cent[ti][0]) * nx + (d.cent[tj][1] - d.cent[ti][1]) * ny < 0.0) {
            nx = -nx; ny = -ny;
        }
        const double de = (d.circ[tj][0] - d.circ[ti][0]) * nx + (d.circ[tj][1] - d.circ[ti][1]) * ny;
        const double dRef = 0.5 * (hi + hj);
        if (de < 0.0) ++d.nonDelaunay;
        if (de < DualMBO::kMinDualDistance * dRef) { ++d.belowFloor; d.clampedEdge[e] = 1; }
        d.dualRatios.push_back(de / dRef);
        d.maxAreaRatio = std::max(d.maxAreaRatio,
                                  std::max(d.area[ti], d.area[tj])
                                      / std::max(1e-300, std::min(d.area[ti], d.area[tj])));
    }
    d.hMean = m.edges.empty() ? 1.0 : hSum / static_cast<double>(m.edges.size());
    std::vector<double> sorted = d.dualRatios;
    std::sort(sorted.begin(), sorted.end());
    d.medianDualRatio = sorted.empty() ? 0.0 : sorted[sorted.size() / 2];

    for (int t = 0; t < NT; ++t) {
        bool ok = true, cl = false;
        for (int k = 0; k < 3; ++k) {
            const int e = m.triangleEdges[t][k];
            if (e < 0 || m.isBoundaryEdge[e]) { ok = false; break; }
            if (d.clampedEdge[e]) cl = true;
        }
        if (ok) { d.tested.push_back(t); d.clamped.push_back(cl ? 1 : 0); }
    }
    for (char c : d.clamped) d.clampedTriangles += c;
    return d;
}

// How far the assembled operator is from annihilating an affine field sampled
// at the cell centres: rho_t = |(M^-1 K phi)_t| * h_mean / gamma, over four
// directions of the gradient, as an rms and a max over the triangles with three
// interior edges. Prop. 4 says the two-point weight gives rounding.
struct Residual { double rms = 0.0, mx = 0.0, rmsAway = 0.0; bool ok = false; };

Residual consistencyResidual(const std::shared_ptr<Mesh> &m, const DualStats &d,
                             DualMBO::PenaltyWeight w, double gammaEval) {
    Residual out;
    const int NT = static_cast<int>(m->triangles.size());
    // No interface alignment and no disk pin: the operator under test is the
    // interior coupling alone, and every triangle read below has three interior
    // edges, so no eliminated row reaches it.
    DualMBO s(m, 1, gammaEval);
    s.setPinDiskCenters(false);
    s.setPenaltyWeight(w);
    s.initialize();
    const Eigen::SparseMatrix<std::complex<double>> &K = s.stiffnessMatrix();
    const Eigen::SparseMatrix<std::complex<double>> &M = s.massMatrix();

    double ss = 0.0, ssAway = 0.0;
    int cnt = 0, cntAway = 0;
    for (int k = 0; k < 4; ++k) {
        const double th = k * M_PI / 4.0;
        const double ax = std::cos(th), ay = std::sin(th);
        Eigen::VectorXcd phi(NT);
        for (int t = 0; t < NT; ++t)
            phi[t] = std::complex<double>(ax * d.circ[t][0] + ay * d.circ[t][1], 0.0);
        const Eigen::VectorXcd r = K * phi;
        for (std::size_t i = 0; i < d.tested.size(); ++i) {
            const int t = d.tested[i];
            const double mtt = M.coeff(t, t).real();
            if (mtt <= 0.0) continue;
            const double rho = std::abs(r[t]) / mtt * d.hMean / gammaEval;
            ss += rho * rho;
            out.mx = std::max(out.mx, rho);
            ++cnt;
            if (!d.clamped[i]) { ssAway += rho * rho; ++cntAway; }
        }
    }
    out.rms = cnt ? std::sqrt(ss / cnt) : 0.0;
    out.rmsAway = cntAway ? std::sqrt(ssAway / cntAway) : 0.0;
    out.ok = true;
    return out;
}

// The disk's answer, known before the run: four +1/4 cones on the diagonals,
// the rotation fixed by the centre pin. Returns the worst deviation of a cone's
// angle from 45 + k*90 degrees, and the largest cone radius.
struct ConeError { double angleDeg = 0.0, radius = 0.0; };

ConeError diskConeError(const Mesh &m, const metrics::SingularityReport &sr) {
    ConeError e;
    for (const metrics::Singularity &sg : sr.interior) {
        const Point d = m.vertices[sg.vertex];
        double ang = std::atan2(d[1], d[0]) * 180.0 / M_PI;
        if (ang < 0.0) ang += 360.0;
        e.angleDeg = std::max(e.angleDeg, std::fabs(std::fmod(ang, 90.0) - 45.0));
        e.radius = std::max(e.radius, normP(d));
    }
    return e;
}

} // namespace

int main(int argc, char **argv) {
    std::string outDir = "results";
    int levels = 4;
    double h = 0.03;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--levels" && i + 1 < argc) levels = std::stoi(argv[++i]);
        else if (a == "--h" && i + 1 < argc) h = std::stod(argv[++i]);
        else if (a == "--weight" && i + 1 < argc) setDefaultPenaltyWeight(parsePenaltyWeight(argv[++i]));
        else if (a == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--out DIR] [--levels N] [--h SIZE] [--weight min|harm|orth]\n";
            return 0;
        }
    }
    ensureDir(outDir);
    Verdicts v;

    banner("E1  Verification and parameter study");
    std::cout << "Sweeps run at h = " << h << "; refinement runs " << levels << " levels.\n"
              << "Energies are the face functional of Sec. 5 at the fixed evaluation "
              << "weight gamma = " << kGammaEval << ".\n"
              << "DualMBO penalty weight for (a)-(d): "
              << penaltyWeightName(defaultPenaltyWeight()) << ".\n";

    // -----------------------------------------------------------------------
    // The four sweep domains, at one size.
    // -----------------------------------------------------------------------
    struct Case { std::string name; std::shared_ptr<Mesh> mesh; bool knownIndex; int expect; };
    std::vector<Case> cases;
    cases.push_back({"square",  squareDomain(h),            true, 0});
    cases.push_back({"disk",    diskDomain(h),              true, 4});
    cases.push_back({"lshape",  lShapeDomain(h),            true, 0});
    cases.push_back({"annulus", annulusDomain(0.4, 1.0, h), true, 0});

    heading("Domains");
    {
        Table t({"domain", "vertices", "triangles", "chi", "boundary quarters"});
        for (const Case &c : cases) {
            const std::vector<int> q = metrics::cornerQuarters(*c.mesh);
            int qs = 0;
            for (std::size_t i = 0; i < q.size(); ++i) if (c.mesh->isBoundaryVertex[i]) qs += q[i];
            t.row({c.name, num((int)c.mesh->vertices.size()), num((int)c.mesh->triangles.size()),
                   num(metrics::eulerCharacteristic(*c.mesh)), num(qs)});
        }
        t.print();
    }

    // -----------------------------------------------------------------------
    // (a) gamma
    // -----------------------------------------------------------------------
    heading("E1(a)  gamma sweep");
    {
        Csv csv(outDir + "/E1_gamma.csv",
                {"domain", "gamma", "triangles", "energy", "singularities", "index4_sum",
                 "ph_residual", "iterations", "converged", "seconds"});
        const std::vector<double> gammas{1.0, 5.0, 10.0, 20.0, 50.0, 100.0};

        for (const Case &c : cases) {
            Table t({"gamma", "energy", "#sing", "sum 4I", "PH", "iters", "s"});
            std::vector<int> counts;
            for (double g : gammas) {
                MethodOptions o;
                o.gamma = g;
                const Solved s = solve(c.mesh, o);
                if (!s.run.ok) { v.warn(c.name + " gamma=" + num(g) + ": " + s.run.error); continue; }
                counts.push_back(static_cast<int>(s.sing.interior.size()));
                t.row({num(g, 4), num(s.energy, 6), num((int)s.sing.interior.size()),
                       num(s.sing.interiorIndex4Sum), num(s.sing.poincareHopfResidual),
                       num(s.run.iterations), num(s.run.totalSeconds, 3)});
                csv.row({{"domain", c.name}, {"gamma", num(g, 4)},
                         {"triangles", num((int)c.mesh->triangles.size())},
                         {"energy", num(s.energy, 10)},
                         {"singularities", num((int)s.sing.interior.size())},
                         {"index4_sum", num(s.sing.interiorIndex4Sum)},
                         {"ph_residual", num(s.sing.poincareHopfResidual)},
                         {"iterations", num(s.run.iterations)},
                         {"converged", num(s.run.converged)},
                         {"seconds", num(s.run.totalSeconds, 4)}});
            }
            std::cout << "\n  " << c.name << "\n";
            t.print(std::cout, "    ");
            const bool stable = !counts.empty() &&
                std::equal(counts.begin() + 1, counts.end(), counts.begin());
            v.check(stable, c.name + ": the singularity count is the same at every gamma");
        }
    }

    // -----------------------------------------------------------------------
    // (b) tau
    // -----------------------------------------------------------------------
    heading("E1(b)  tau sweep about the D^2/10 heuristic");
    {
        Csv csv(outDir + "/E1_tau.csv",
                {"domain", "tau_scale", "energy", "singularities", "index4_sum",
                 "ph_residual", "iterations", "converged"});
        const std::vector<double> scales{0.1, 0.5, 1.0, 2.0, 10.0};

        for (const Case &c : cases) {
            Table t({"tau/(D^2/10)", "energy", "#sing", "sum 4I", "PH", "iters"});
            std::vector<int> counts;
            for (double s0 : scales) {
                MethodOptions o;
                o.tauScale = s0;
                const Solved s = solve(c.mesh, o);
                if (!s.run.ok) { v.warn(c.name + " tau x" + num(s0) + ": " + s.run.error); continue; }
                counts.push_back(static_cast<int>(s.sing.interior.size()));
                t.row({num(s0, 3), num(s.energy, 6), num((int)s.sing.interior.size()),
                       num(s.sing.interiorIndex4Sum), num(s.sing.poincareHopfResidual),
                       num(s.run.iterations)});
                csv.row({{"domain", c.name}, {"tau_scale", num(s0, 4)},
                         {"energy", num(s.energy, 10)},
                         {"singularities", num((int)s.sing.interior.size())},
                         {"index4_sum", num(s.sing.interiorIndex4Sum)},
                         {"ph_residual", num(s.sing.poincareHopfResidual)},
                         {"iterations", num(s.run.iterations)},
                         {"converged", num(s.run.converged)}});
            }
            std::cout << "\n  " << c.name << "\n";
            t.print(std::cout, "    ");
            const bool stable = !counts.empty() &&
                std::equal(counts.begin() + 1, counts.end(), counts.begin());
            if (!stable) v.warn(c.name + ": the singularity count depends on tau over x0.1..x10");
            else std::cout << "    (count stable over the whole sweep)\n";
        }
    }

    // -----------------------------------------------------------------------
    // (a/b) the product identity
    // -----------------------------------------------------------------------
    heading("E1(a/b)  gamma and tau enter only through their product");
    {
        Csv csv(outDir + "/E1_product.csv",
                {"domain", "gamma_a", "tau_a", "gamma_b", "tau_b", "product",
                 "max_abs_difference", "energy_a", "energy_b"});
        // Four (gamma, tau) pairs at the same product, spread over two decades
        // of each. If the identity holds, all four give the same field.
        const std::vector<std::pair<double, double>> pairs{
            {10.0, 1.0}, {1.0, 10.0}, {100.0, 0.1}, {2.0, 5.0}};

        for (const Case &c : cases) {
            Solved ref;
            double worst = 0.0;
            for (std::size_t k = 0; k < pairs.size(); ++k) {
                MethodOptions o;
                o.gamma = pairs[k].first;
                o.tauScale = pairs[k].second;
                const Solved s = solve(c.mesh, o);
                if (!s.run.ok) { v.warn(c.name + ": product check solve failed"); break; }
                if (k == 0) { ref = s; continue; }
                double d = 0.0;
                for (int i = 0; i < s.run.u.size(); ++i)
                    d = std::max(d, std::abs(s.run.u[i] - ref.run.u[i]));
                worst = std::max(worst, d);
                csv.row({{"domain", c.name}, {"gamma_a", num(pairs[0].first, 4)},
                         {"tau_a", num(pairs[0].second, 4)},
                         {"gamma_b", num(pairs[k].first, 4)},
                         {"tau_b", num(pairs[k].second, 4)},
                         {"product", num(pairs[k].first * pairs[k].second, 6)},
                         {"max_abs_difference", num(d, 6)},
                         {"energy_a", num(ref.energy, 10)},
                         {"energy_b", num(s.energy, 10)}});
            }
            std::cout << "  " << c.name << ": worst |u(gamma,tau) - u(10,1)| over three "
                      << "equal-product pairs = " << num(worst, 3) << "\n";
            v.check(worst < 1e-9,
                    c.name + ": the field depends on gamma and tau only through gamma*tau");
        }
    }

    // -----------------------------------------------------------------------
    // (c) h-refinement, on the disk
    // -----------------------------------------------------------------------
    heading("E1(c)  h-refinement on the disk");
    {
        Csv csv(outDir + "/E1_refine.csv",
                {"domain", "h", "vertices", "triangles", "energy", "singularities",
                 "index4_sum", "ph_residual", "iterations", "seconds",
                 "max_cone_angle_error_deg", "max_cone_radius_frac"});
        Table t({"h", "triangles", "energy", "#sing", "sum 4I", "PH",
                 "max |angle mod 90 - 45| deg", "max r/R", "s"});

        std::vector<double> energies, hs;
        for (double hi : refinementSizes(levels)) {
            std::shared_ptr<Mesh> m = diskDomain(hi);
            MethodOptions o;
            const Solved s = solve(m, o);
            if (!s.run.ok) { v.warn("disk h=" + num(hi) + ": " + s.run.error); continue; }

            // The disk's answer, known before the run: four +1/4 cones on the
            // diagonals. The centre pin is what fixes the rotation -- without it
            // the angles are arbitrary and this column means nothing.
            const ConeError ce = diskConeError(*m, s.sing);
            const double worstAngle = ce.angleDeg, worstRadius = ce.radius;

            hs.push_back(hi);
            energies.push_back(s.energy);
            t.row({num(hi, 4), num((int)m->triangles.size()), num(s.energy, 6),
                   num((int)s.sing.interior.size()), num(s.sing.interiorIndex4Sum),
                   num(s.sing.poincareHopfResidual), num(worstAngle, 4), num(worstRadius, 4),
                   num(s.run.totalSeconds, 3)});
            csv.row({{"domain", "disk"}, {"h", num(hi, 6)},
                     {"vertices", num((int)m->vertices.size())},
                     {"triangles", num((int)m->triangles.size())},
                     {"energy", num(s.energy, 10)},
                     {"singularities", num((int)s.sing.interior.size())},
                     {"index4_sum", num(s.sing.interiorIndex4Sum)},
                     {"ph_residual", num(s.sing.poincareHopfResidual)},
                     {"iterations", num(s.run.iterations)},
                     {"seconds", num(s.run.totalSeconds, 4)},
                     {"max_cone_angle_error_deg", num(worstAngle, 6)},
                     {"max_cone_radius_frac", num(worstRadius, 6)}});

            v.check(s.sing.poincareHopfResidual == 0,
                    "disk h=" + num(hi, 4) + ": Poincare-Hopf holds");
            v.check(s.sing.interiorIndex4Sum == 4,
                    "disk h=" + num(hi, 4) + ": the interior index total is +1");
        }
        t.print();

        // The energy of a cross field with singularities diverges as h -> 0 --
        // the Ginzburg-Landau |log eps| -- so what converges is not the energy
        // but the *rate*. Fit E = a + s log(1/h) by least squares and report s;
        // the singularity-free domains below are the control, where E stays at
        // rounding whatever h is.
        if (hs.size() >= 2) {
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            const double n = static_cast<double>(hs.size());
            for (std::size_t i = 0; i < hs.size(); ++i) {
                const double x = std::log(1.0 / hs[i]), y = energies[i];
                sx += x; sy += y; sxx += x * x; sxy += x * y;
            }
            const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
            std::cout << "    E = a + s log(1/h) fitted over " << hs.size()
                      << " levels: s = " << num(slope, 5) << "\n";
            v.check(slope > 0.0,
                    "disk: the energy grows like log(1/h), as four cones require");
        }

        // The same refinement on the square, the L and the annulus, which is
        // where the index total being h-independent is the interesting claim.
        for (const Case &c : cases) {
            if (c.name == "disk") continue;
            for (double hi : refinementSizes(std::min(levels, 3))) {
                std::shared_ptr<Mesh> m =
                    c.name == "square"  ? squareDomain(hi)
                  : c.name == "lshape"  ? lShapeDomain(hi)
                                        : annulusDomain(0.4, 1.0, hi);
                MethodOptions o;
                const Solved s = solve(m, o);
                if (!s.run.ok) continue;
                csv.row({{"domain", c.name}, {"h", num(hi, 6)},
                         {"vertices", num((int)m->vertices.size())},
                         {"triangles", num((int)m->triangles.size())},
                         {"energy", num(s.energy, 10)},
                         {"singularities", num((int)s.sing.interior.size())},
                         {"index4_sum", num(s.sing.interiorIndex4Sum)},
                         {"ph_residual", num(s.sing.poincareHopfResidual)},
                         {"iterations", num(s.run.iterations)},
                         {"seconds", num(s.run.totalSeconds, 4)},
                         {"max_cone_angle_error_deg", ""},
                         {"max_cone_radius_frac", ""}});
                v.check(s.sing.poincareHopfResidual == 0,
                        c.name + " h=" + num(hi, 4) + ": Poincare-Hopf holds");
                if (c.knownIndex)
                    v.check(s.sing.interiorIndex4Sum == c.expect,
                            c.name + " h=" + num(hi, 4) + ": interior index total is the expected "
                            + num(c.expect));
            }
        }
    }

    // -----------------------------------------------------------------------
    // (d) MBO convergence
    // -----------------------------------------------------------------------
    heading("E1(d)  MBO convergence and the Lyapunov energy");
    {
        Csv csv(outDir + "/E1_convergence.csv",
                {"domain", "iteration", "increment", "energy", "dualmbo_energy"});
        Table t({"domain", "steps", "increment first -> last", "E_common first -> last",
                 "worst rise (rel)", "E_DualMBO first -> last", "worst rise (rel)"});

        for (const Case &c : cases) {
            MethodOptions o;
            o.recordHistory = true;
            o.forceSteps = true;      // a fixed budget, so the curve has a shape
            o.maxSteps = 40;
            // The single-tau convergence history, for the same reason `solve`
            // forces it: a continuation restarts the increment at every level,
            // so its history is a sawtooth of restarts and says nothing about
            // whether one MBO iteration converges.
            o.tauContinuation = false;
            const FieldRun r = runDualMBO(c.mesh, o);
            if (!r.ok) { v.warn(c.name + ": " + r.error); continue; }
            if (r.energyHistory.empty()) continue;

            for (std::size_t k = 0; k < r.energyHistory.size(); ++k) {
                csv.row({{"domain", c.name}, {"iteration", num((int)k + 1)},
                         {"increment", num(r.incrementHistory[k], 10)},
                         {"energy", num(r.energyHistory[k], 10)},
                         {"dualmbo_energy", num(r.dualMBOEnergyHistory[k], 10)}});
            }

            // The three monotonicity questions, answered rather than assumed.
            auto worstRise = [](const std::vector<double> &e) {
                double w = 0.0;
                for (std::size_t k = 1; k < e.size(); ++k) w = std::max(w, e[k] - e[k - 1]);
                return w;
            };
            const double scaleC = std::max(1e-30, std::fabs(r.energyHistory.front()));
            const double scaleS = std::max(1e-30, std::fabs(r.dualMBOEnergyHistory.front()));
            const double riseC = worstRise(r.energyHistory) / scaleC;
            const double riseS = worstRise(r.dualMBOEnergyHistory) / scaleS;

            t.row({c.name, num((int)r.energyHistory.size()),
                   num(r.incrementHistory.front(), 3) + " -> " + num(r.incrementHistory.back(), 3),
                   num(r.energyHistory.front(), 6) + " -> " + num(r.energyHistory.back(), 6),
                   num(riseC, 3),
                   num(r.dualMBOEnergyHistory.front(), 6) + " -> " + num(r.dualMBOEnergyHistory.back(), 6),
                   num(riseS, 3)});

            // What is asserted is the thing the scheme is actually run on: the
            // increment goes to zero, so the iteration has a fixed point and the
            // stopping criterion means something.
            v.check(r.incrementHistory.back() <= r.incrementHistory.front() + 1e-12,
                    c.name + ": the MBO increment does not grow");
            v.check(r.incrementHistory.back() < 1e-6 * std::max(1.0, r.incrementHistory.front()),
                    c.name + ": the increment has fallen to rounding by step "
                    + num((int)r.incrementHistory.size()));
        }
        t.print();

        // What the shipped stopping criterion costs. The pipeline stops at
        // `error < 2 NT 1e-5`, which on these domains fires after two or three
        // steps -- long before the increment is anywhere near zero -- and the
        // field it stops at is not the fixed point. This is the measurement of
        // that gap, and it is why the field-quality experiments run tighter.
        std::cout << "\n  What the shipped stopping criterion stops at:\n";
        Table st({"domain", "steps at 1e-5", "E_common", "steps at 1e-9", "E_common",
                  "relative gap"});
        for (const Case &c : cases) {
            MethodOptions loose;
            const Solved a = solve(c.mesh, loose);
            MethodOptions tight;
            tight.convergenceTol = 1e-9;
            const Solved b = solve(c.mesh, tight);
            if (!a.run.ok || !b.run.ok) continue;
            const double gap = std::fabs(a.energy - b.energy) / std::max(1e-30, std::fabs(b.energy));
            st.row({c.name, num(a.run.iterations), num(a.energy, 6),
                    num(b.run.iterations), num(b.energy, 6), num(gap, 3)});
        }
        st.print();

        std::cout <<
            "\n  Reading: MBO is monotone in the Esedoglu-Otto heat content, not in a\n"
            "  Dirichlet-type energy, so a rise of a fraction of a percent in either\n"
            "  column above is the scheme behaving correctly. Prop. 2 should be stated\n"
            "  in the heat content and this table is the evidence for which functional\n"
            "  the statement can be made about.\n";
    }

    // -----------------------------------------------------------------------
    // (e) consistency of the p = 0 operator, over the corpus
    // -----------------------------------------------------------------------
    heading("E1(e)  Consistency of the edge weight on an affine field");
    {
        std::cout <<
            "  phi(x) = a.x is sampled at the circumcenters, r = M^-1 K phi is read on\n"
            "  every triangle whose three edges are interior, and rho_t = |r_t| h_mean /\n"
            "  gamma is reported as an rms and a max over four directions of a. A\n"
            "  consistent second-order operator gives rounding. d_e is the signed\n"
            "  distance between the two circumcenters across an edge and d_ref =\n"
            "  (h_i + h_j)/2; d_e/d_ref = 2/3 on an equilateral pair, d_e = 0 on a\n"
            "  cocircular one, d_e < 0 on a non-Delaunay one.\n";

        Csv csv(outDir + "/E1_consistency.csv",
                {"set", "model", "triangles", "interior_edges", "tested_triangles",
                 "non_delaunay_edges", "edges_below_floor", "clamped_triangles",
                 "median_dual_ratio", "max_area_ratio", "weight",
                 "residual_rms", "residual_max", "residual_rms_away_from_floor"});

        const std::vector<DualMBO::PenaltyWeight> weights{
            DualMBO::PenaltyWeight::MinHeight,
            DualMBO::PenaltyWeight::HarmonicHeight,
            DualMBO::PenaltyWeight::Orthogonal};

        for (const std::string set : {"singlemat", "multimat"}) {
            const std::vector<Model> models = corpus(set);
            if (models.empty()) {
                v.warn("E1(e): no models under " + std::string(PAPER_MESH_DIR) + "/" + set);
                continue;
            }
            std::cout << "\n  data/meshes/" << set << "\n";
            Table t({"model", "triangles", "non-Delaunay", "below floor", "median d_e/d_ref",
                     "max area ratio", "rms min", "rms harm", "rms orth", "max orth",
                     "rms orth away from floor"});

            std::map<std::string, double> rmsSum, rmsWorst, rmsAwaySum, rmsAwayWorst;
            int n = 0, orthAtRounding = 0, interiorEdgesTotal = 0, nonDelaunayTotal = 0,
                belowFloorTotal = 0, clampedTrianglesTotal = 0;
            std::vector<double> allRatios;

            for (const Model &mm : models) {
                std::string why;
                std::shared_ptr<Mesh> m = load(mm, why);
                if (!m) { v.warn(mm.name + ": " + why); continue; }
                const int NT = static_cast<int>(m->triangles.size());

                const DualStats d = dualStats(*m);
                allRatios.insert(allRatios.end(), d.dualRatios.begin(), d.dualRatios.end());

                std::map<std::string, Residual> res;
                for (DualMBO::PenaltyWeight w : weights) {
                    const std::string wn = penaltyWeightName(w);
                    Residual r;
                    try {
                        r = consistencyResidual(m, d, w, kGammaEval);
                    } catch (const std::exception &ex) {
                        v.warn(mm.name + "/" + wn + ": " + ex.what());
                        continue;
                    }
                    res[wn] = r;
                    rmsSum[wn] += r.rms;
                    rmsWorst[wn] = std::max(rmsWorst[wn], r.rms);
                    rmsAwaySum[wn] += r.rmsAway;
                    rmsAwayWorst[wn] = std::max(rmsAwayWorst[wn], r.rmsAway);
                    csv.row({{"set", set}, {"model", mm.name}, {"triangles", num(NT)},
                             {"interior_edges", num(d.interiorEdges)},
                             {"tested_triangles", num((int)d.tested.size())},
                             {"non_delaunay_edges", num(d.nonDelaunay)},
                             {"edges_below_floor", num(d.belowFloor)},
                             {"clamped_triangles", num(d.clampedTriangles)},
                             {"median_dual_ratio", num(d.medianDualRatio, 6)},
                             {"max_area_ratio", num(d.maxAreaRatio, 6)},
                             {"weight", wn}, {"residual_rms", num(r.rms, 6)},
                             {"residual_max", num(r.mx, 6)},
                             {"residual_rms_away_from_floor", num(r.rmsAway, 6)}});
                }
                if (res.size() != weights.size()) continue;
                ++n;
                interiorEdgesTotal += d.interiorEdges;
                nonDelaunayTotal += d.nonDelaunay;
                belowFloorTotal += d.belowFloor;
                clampedTrianglesTotal += d.clampedTriangles;
                if (res["orth"].rms < 1e-12) ++orthAtRounding;
                t.row({mm.name, num(NT), num(d.nonDelaunay), num(d.belowFloor),
                       num(d.medianDualRatio, 4),
                       num(d.maxAreaRatio, 3), num(res["min"].rms, 3), num(res["harm"].rms, 3),
                       num(res["orth"].rms, 3), num(res["orth"].mx, 3), num(res["orth"].rmsAway, 3)});
                // Prop. 4 is a statement about the unclamped operator: assert it
                // away from the floor, and let the table show what the clamp costs.
                v.check(res["orth"].rmsAway < 1e-10,
                        mm.name + ": the two-point weight annihilates the affine field away from "
                        + num(d.clampedTriangles) + " floor-clamped triangle(s)");
            }
            t.print();
            if (n) {
                std::sort(allRatios.begin(), allRatios.end());
                std::cout << "  " << n << " model(s), " << interiorEdgesTotal << " interior edges: "
                          << nonDelaunayTotal << " non-Delaunay, " << belowFloorTotal
                          << " below the floor d_e < " << DualMBO::kMinDualDistance
                          << " d_ref (" << clampedTrianglesTotal
                          << " triangles touch one); median d_e/d_ref over all edges "
                          << num(allRatios[allRatios.size() / 2], 4) << "\n";
                std::cout << "  mean rms residual: min " << num(rmsSum["min"] / n, 3)
                          << ", harm " << num(rmsSum["harm"] / n, 3)
                          << ", orth " << num(rmsSum["orth"] / n, 3)
                          << "; worst rms: min " << num(rmsWorst["min"], 3)
                          << ", harm " << num(rmsWorst["harm"], 3)
                          << ", orth " << num(rmsWorst["orth"], 3) << "\n";
                std::cout << "  away from floor-clamped edges: mean rms min "
                          << num(rmsAwaySum["min"] / n, 3) << ", harm " << num(rmsAwaySum["harm"] / n, 3)
                          << ", orth " << num(rmsAwaySum["orth"] / n, 3)
                          << "; worst orth " << num(rmsAwayWorst["orth"], 3) << "\n";
                std::cout << "  models with the two-point residual below 1e-12 on every tested "
                          << "triangle: " << orthAtRounding << "/" << n << "\n";
            }
        }
    }

    // -----------------------------------------------------------------------
    // (f) does the inconsistency matter? -- the graded-mesh sweep
    // -----------------------------------------------------------------------
    heading("E1(f)  Graded meshes: what the inconsistency costs the field");
    {
        std::cout <<
            "  The corpus is near-uniform, so (e)'s O(1) residual can be a defect that is\n"
            "  provable and invisible. This grades the disk -- concentric rings whose radial\n"
            "  spacing grows by `growth` per ring, constrained Delaunay, no Steiner points --\n"
            "  and reads two things at each grading: the operator's residual on an affine\n"
            "  field, and where the four cones actually land. The disk is the domain whose\n"
            "  answer is known in closed form (four +1/4 cones at 45 + k*90 degrees from the\n"
            "  centre, fixed by the centre pin), so the cone column is *error*, not a proxy.\n";

        Csv csv(outDir + "/E1_graded.csv",
                {"growth", "vertices", "triangles", "max_area_ratio", "median_dual_ratio",
                 "non_delaunay_edges", "edges_below_floor", "weight",
                 "residual_rms", "residual_max", "residual_rms_away_from_floor",
                 "singularities", "index4_sum", "ph_residual",
                 "max_cone_angle_error_deg", "max_cone_radius_frac",
                 "energy_twopoint", "energy_minheight", "iterations"});

        const std::vector<double> growths{1.0, 1.05, 1.10, 1.15, 1.20};
        const std::vector<DualMBO::PenaltyWeight> weights{
            DualMBO::PenaltyWeight::MinHeight,
            DualMBO::PenaltyWeight::Orthogonal};

        Table t({"growth", "triangles", "max area ratio", "median d_e/d_ref", "clamped edges",
                 "weight", "residual rms", "rms away from floor", "#sing", "sum 4I", "PH",
                 "max cone angle err (deg)", "max r/R", "E (two-point)"});
        std::map<std::string, std::vector<double>> coneErr;

        for (double g : growths) {
            std::shared_ptr<Mesh> m;
            try {
                m = gradedDiskDomain(1.0, 50, g);
            } catch (const std::exception &ex) {
                v.warn("graded disk growth=" + num(g, 3) + ": " + ex.what());
                continue;
            }
            const DualStats d = dualStats(*m);
            const metrics::EdgeWeights wTP = metrics::edgeWeightsTwoPoint(*m, kGammaEval);
            const metrics::EdgeWeights wMin = metrics::edgeWeights(*m, kGammaEval);

            for (DualMBO::PenaltyWeight w : weights) {
                const std::string wn = penaltyWeightName(w);
                Residual r;
                try {
                    r = consistencyResidual(m, d, w, kGammaEval);
                } catch (const std::exception &ex) {
                    v.warn("graded disk growth=" + num(g, 3) + "/" + wn + ": " + ex.what());
                    continue;
                }

                MethodOptions o;
                o.penaltyWeight = w;
                o.convergenceTol = 1e-9;
                const Solved s = solve(m, o);
                if (!s.run.ok) {
                    v.warn("graded disk growth=" + num(g, 3) + "/" + wn + ": " + s.run.error);
                    continue;
                }
                const ConeError ce = diskConeError(*m, s.sing);
                const double eTP = metrics::commonEnergy(*m, wTP, s.run.u);
                const double eMin = metrics::commonEnergy(*m, wMin, s.run.u);
                coneErr[wn].push_back(ce.angleDeg);

                t.row({num(g, 3), num((int)m->triangles.size()), num(d.maxAreaRatio, 3),
                       num(d.medianDualRatio, 3), num(d.belowFloor), wn, num(r.rms, 3),
                       num(r.rmsAway, 3),
                       num((int)s.sing.interior.size()), num(s.sing.interiorIndex4Sum),
                       num(s.sing.poincareHopfResidual), num(ce.angleDeg, 4),
                       num(ce.radius, 4), num(eTP, 6)});
                csv.row({{"growth", num(g, 4)},
                         {"vertices", num((int)m->vertices.size())},
                         {"triangles", num((int)m->triangles.size())},
                         {"max_area_ratio", num(d.maxAreaRatio, 6)},
                         {"median_dual_ratio", num(d.medianDualRatio, 6)},
                         {"non_delaunay_edges", num(d.nonDelaunay)},
                         {"edges_below_floor", num(d.belowFloor)},
                         {"weight", wn},
                         {"residual_rms", num(r.rms, 6)}, {"residual_max", num(r.mx, 6)},
                         {"residual_rms_away_from_floor", num(r.rmsAway, 6)},
                         {"singularities", num((int)s.sing.interior.size())},
                         {"index4_sum", num(s.sing.interiorIndex4Sum)},
                         {"ph_residual", num(s.sing.poincareHopfResidual)},
                         {"max_cone_angle_error_deg", num(ce.angleDeg, 6)},
                         {"max_cone_radius_frac", num(ce.radius, 6)},
                         {"energy_twopoint", num(eTP, 10)},
                         {"energy_minheight", num(eMin, 10)},
                         {"iterations", num(s.run.iterations)}});

                v.check(s.sing.poincareHopfResidual == 0,
                        "graded disk growth=" + num(g, 3) + "/" + wn + ": Poincare-Hopf holds");
            }
        }
        t.print();

        // The one comparison the sweep exists to make, stated rather than left
        // to the reader: does the cone placement error grow with the grading for
        // the inconsistent weight and not for the consistent one?
        // The comparison is *paired*: both weights run on the same mesh at each
        // level, so a difference between them is the operator and nothing else.
        {
            const std::vector<double> &a = coneErr["min"];
            const std::vector<double> &b = coneErr["orth"];
            int betterOrth = 0, levels = 0;
            double worstA = 0.0, worstB = 0.0;
            for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) {
                ++levels;
                if (b[i] < a[i]) ++betterOrth;
                worstA = std::max(worstA, a[i]);
                worstB = std::max(worstB, b[i]);
            }
            std::cout << "\n  Paired over " << levels << " grading level(s), on identical meshes:\n"
                      << "    the two-point weight places the cones closer to their known angles on "
                      << betterOrth << "/" << levels << ",\n"
                      << "    worst cone angle error min-height " << num(worstA, 4)
                      << " deg against two-point " << num(worstB, 4) << " deg.\n"
                      << "\n  This is the downstream half of Prop. 4. The corpus is near-uniform and\n"
                      << "  there the two weights give almost the same field; grade the mesh and the\n"
                      << "  inconsistency stops being invisible. Where the two columns agree the\n"
                      << "  honest reading is that the defect is real at the operator level and does\n"
                      << "  not reach the field on that mesh.\n";
            v.check(betterOrth >= levels - 1,
                    "graded disk: the consistent weight places the cones at least as well at every"
                    " grading level but at most one");
        }
    }

    banner("E1 summary");
    std::cout << (v.failures ? kFail : kPass) << " " << v.failures << " failure(s), "
              << v.warnings << " warning(s). CSVs in " << outDir << "/\n";
    return v.failures ? 1 : 0;
}
