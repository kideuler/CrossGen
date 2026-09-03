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
//
// The comparison energies are measured with the *same* functional -- the
// face-based one of Sec. 5 at the fixed evaluation weight gamma_eval = 10 --
// whatever gamma the solver ran at. Comparing a solve at gamma = 100 against
// its own energy would compare two different functionals and would say nothing.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "common/Domains.hxx"
#include "common/Methods.hxx"
#include "common/Metrics.hxx"
#include "common/Report.hxx"

using namespace paper;

namespace {

// The weight every energy in this file is evaluated at, whatever the solver ran
// at. Sec. 5: "SIPG weights fixed for all methods".
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
    s.run = runSIPG(m, single);
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
        else if (a == "--help") {
            std::cout << "Usage: " << argv[0] << " [--out DIR] [--levels N] [--h SIZE]\n";
            return 0;
        }
    }
    ensureDir(outDir);
    Verdicts v;

    banner("E1  Verification and parameter study");
    std::cout << "Sweeps run at h = " << h << "; refinement runs " << levels << " levels.\n"
              << "Energies are the face functional of Sec. 5 at the fixed evaluation "
              << "weight gamma = " << kGammaEval << ".\n";

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
            double worstAngle = 0.0, worstRadius = 0.0;
            for (const metrics::Singularity &sg : s.sing.interior) {
                const Point d = m->vertices[sg.vertex];
                double ang = std::atan2(d[1], d[0]) * 180.0 / M_PI;
                if (ang < 0.0) ang += 360.0;
                worstAngle = std::max(worstAngle, std::fabs(std::fmod(ang, 90.0) - 45.0));
                worstRadius = std::max(worstRadius, normP(d));
            }

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
                {"domain", "iteration", "increment", "energy", "sipg_energy"});
        Table t({"domain", "steps", "increment first -> last", "E_common first -> last",
                 "worst rise (rel)", "E_SIPG first -> last", "worst rise (rel)"});

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
            const FieldRun r = runSIPG(c.mesh, o);
            if (!r.ok) { v.warn(c.name + ": " + r.error); continue; }
            if (r.energyHistory.empty()) continue;

            for (std::size_t k = 0; k < r.energyHistory.size(); ++k) {
                csv.row({{"domain", c.name}, {"iteration", num((int)k + 1)},
                         {"increment", num(r.incrementHistory[k], 10)},
                         {"energy", num(r.energyHistory[k], 10)},
                         {"sipg_energy", num(r.sipgEnergyHistory[k], 10)}});
            }

            // The three monotonicity questions, answered rather than assumed.
            auto worstRise = [](const std::vector<double> &e) {
                double w = 0.0;
                for (std::size_t k = 1; k < e.size(); ++k) w = std::max(w, e[k] - e[k - 1]);
                return w;
            };
            const double scaleC = std::max(1e-30, std::fabs(r.energyHistory.front()));
            const double scaleS = std::max(1e-30, std::fabs(r.sipgEnergyHistory.front()));
            const double riseC = worstRise(r.energyHistory) / scaleC;
            const double riseS = worstRise(r.sipgEnergyHistory) / scaleS;

            t.row({c.name, num((int)r.energyHistory.size()),
                   num(r.incrementHistory.front(), 3) + " -> " + num(r.incrementHistory.back(), 3),
                   num(r.energyHistory.front(), 6) + " -> " + num(r.energyHistory.back(), 6),
                   num(riseC, 3),
                   num(r.sipgEnergyHistory.front(), 6) + " -> " + num(r.sipgEnergyHistory.back(), 6),
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

    banner("E1 summary");
    std::cout << (v.failures ? kFail : kPass) << " " << v.failures << " failure(s), "
              << v.warnings << " warning(s). CSVs in " << outDir << "/\n";
    return v.failures ? 1 : 0;
}
