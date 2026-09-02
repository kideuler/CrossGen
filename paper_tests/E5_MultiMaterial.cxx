// E5 -- Multi-material domains (Sec. 6, 0.65 page). Tier 1.
//
// This experiment carries claim 4 and half the paper's motivation, so it is in
// three parts and the last two are the ones with numbers in them.
//
//   Part 1  The field at the interfaces. Interface alignment error per method,
//           over every interface edge, taken one side at a time -- because the
//           two sides of an interface carry different crosses by design and an
//           average of them would measure neither.
//
//   Part 2  The junction vertices, which is where the claim actually lives.
//           At a vertex where k >= 3 materials meet, a vertex-based field has
//           one value and k tangency constraints, and no cross satisfies more
//           than two of them unless the sectors happen to be right angles. The
//           **junction residual** is the number that says so: for a vertex
//           field, the largest angle between the vertex's one cross and the
//           tangent of an incident interface edge; for a face field, the
//           largest over incident faces of the angle between *that face's*
//           cross and the interfaces bounding it. Face DOFs make it a
//           per-face question and it goes to zero; vertex DOFs make it one
//           question with k answers and it does not.
//
//   Part 3  The layout. Mixed-zone count (R1, and the number the application
//           audience cares about most: it should be 0), irregular vertices per
//           material region rather than globally (R4 -- a badly placed one
//           inside a thin layer is worse than one in the bulk), and the
//           quality and CFL-proxy columns E4 uses.
//
//   Rider   E5(a), the boundary-condition ablation: hard-pinned Dirichlet rows
//           against the weak Nitsche form the same assembly already builds.
//           Sec. 4.4 asks for this to be presented as a choice with a
//           measurement behind it, and for the results to say honestly which
//           one they use.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "common/Corpus.hxx"
#include "common/Domains.hxx"
#include "common/Layout.hxx"
#include "common/Methods.hxx"
#include "common/Metrics.hxx"
#include "common/Report.hxx"

using namespace paper;

namespace {

constexpr double kGammaEval = 10.0;
constexpr double kTol = 1e-9;

std::vector<std::string> split(const std::string &s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, sep)) if (!tok.empty()) out.push_back(tok);
    return out;
}
double mean(const std::vector<double> &v) {
    return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
double deg(double rad) { return rad * 180.0 / M_PI; }

// Vertices where three or more materials meet, or where an interface ends on
// dS. These are the vertices Sec. 4.4's proposition is about.
struct Junctions {
    std::vector<int> vertices;
    std::vector<int> materialsAt;    // parallel: how many distinct materials
    std::vector<std::vector<int>> edgesAt;  // the interface edges at each
    // How far the junction's sectors are from multiples of a right angle, in
    // degrees, worst sector. This is the number that decides whether the
    // junction is a genuine conflict or a control: a cross is invariant under a
    // quarter turn, so a junction whose sectors are all 90 or 180 degrees can be
    // satisfied by *one* cross and a vertex-based field has no difficulty there.
    // Only an oblique junction forces the k constraints apart, and then it
    // forces them apart by about half the obliquity.
    std::vector<double> obliquity;
};

Junctions junctions(const Mesh &m) {
    Junctions j;
    const std::vector<int> ie = interfaceEdges(m);
    std::map<int, std::vector<int>> at;
    for (int e : ie) {
        at[m.edges[e][0]].push_back(e);
        at[m.edges[e][1]].push_back(e);
    }
    for (const auto &[v, edges] : at) {
        std::set<int> mats;
        const auto &vt = m.vertexTriangles;
        for (int k = vt.rowPtr[v]; k < vt.rowPtr[v + 1]; ++k)
            mats.insert(m.triangleMatId[vt.colIdx[k]]);
        // A junction: three or more interface edges meeting, three or more
        // materials incident, or an interface landing on dS (where the
        // interface tangent and the boundary tangent are two constraints on
        // whatever sits in the corner).
        const bool onBoundary = m.isBoundaryVertex[v];
        if (edges.size() >= 3 || mats.size() >= 3 || (onBoundary && !edges.empty())) {
            j.vertices.push_back(v);
            j.materialsAt.push_back(static_cast<int>(mats.size()));
            j.edgesAt.push_back(edges);

            // The sectors: the rays of every constraint at v -- the interface
            // edges, and the boundary edges too where v is on dS, since dS is a
            // tangency constraint in exactly the same way.
            std::vector<int> rays = edges;
            if (onBoundary) {
                const auto &vtx = m.vertexTriangles;
                for (int k = vtx.rowPtr[v]; k < vtx.rowPtr[v + 1]; ++k) {
                    for (int c = 0; c < 3; ++c) {
                        const int e = m.triangleEdges[vtx.colIdx[k]][c];
                        if (!m.isBoundaryEdge[e]) continue;
                        if (m.edges[e][0] == v || m.edges[e][1] == v) rays.push_back(e);
                    }
                }
            }
            std::vector<double> ang;
            for (int e : rays) {
                const int w = (m.edges[e][0] == v) ? m.edges[e][1] : m.edges[e][0];
                const Point d = m.vertices[w] - m.vertices[v];
                if (normP(d) > 1e-14) ang.push_back(std::atan2(d[1], d[0]));
            }
            std::sort(ang.begin(), ang.end());
            ang.erase(std::unique(ang.begin(), ang.end()), ang.end());
            double worst = 0.0;
            for (std::size_t i = 0; i < ang.size(); ++i) {
                double sec = ang[(i + 1) % ang.size()] - ang[i];
                if (sec <= 0.0) sec += 2.0 * M_PI;
                const double q = sec * 2.0 / M_PI;       // in quarter turns
                worst = std::max(worst, std::fabs(q - std::round(q)) * 90.0);
            }
            j.obliquity.push_back(ang.size() < 2 ? 0.0 : worst);
        }
    }
    return j;
}

// The junction residual, for a face field: the worst, over the junction's
// interface edges and the faces on either side of each, of the angle between
// that face's cross and the edge.
double junctionResidualFace(const Mesh &m, const Eigen::VectorXcd &u,
                            const std::vector<int> &edges) {
    double worst = 0.0;
    for (int e : edges) {
        const Point d = m.vertices[m.edges[e][1]] - m.vertices[m.edges[e][0]];
        if (normP(d) < 1e-14) continue;
        const double theta = std::atan2(d[1], d[0]);
        for (int s = 0; s < 2; ++s) {
            const int t = m.edgeTriangles[e][s];
            if (t < 0) continue;
            worst = std::max(worst, metrics::crossToDirection(u[t], theta));
        }
    }
    return worst;
}

// The same for a vertex field: one cross at the junction, every incident
// interface edge asking it to be tangent. This is the quantitative form of "a
// continuous vertex-based field is ill posed at a junction".
double junctionResidualVertex(const Mesh &m, const Eigen::VectorXcd &uv, int v,
                              const std::vector<int> &edges) {
    if (v >= uv.size()) return 0.0;
    double worst = 0.0;
    for (int e : edges) {
        const Point d = m.vertices[m.edges[e][1]] - m.vertices[m.edges[e][0]];
        if (normP(d) < 1e-14) continue;
        worst = std::max(worst, metrics::crossToDirection(uv[v], std::atan2(d[1], d[0])));
    }
    return worst;
}

// Irregular vertices per material, printed compactly.
std::string perMaterial(const std::vector<int> &v) {
    if (v.empty()) return "-";
    std::ostringstream oss;
    for (std::size_t i = 0; i < v.size(); ++i) oss << (i ? "/" : "") << v[i];
    return oss.str();
}

struct Case {
    std::string name;
    std::shared_ptr<Mesh> mesh;
    std::string note;
};

} // namespace

int main(int argc, char **argv) {
    std::string outDir = "results";
    std::string objDir, vtkDir;
    std::vector<std::string> methods{"SIPG", "B1"};
    std::vector<std::string> only;
    double target = 0.05;
    double h = 0.04;
    bool doLayout = true;
    bool doJunctionDomains = true;
    bool diskTemplates = false;
    int limit = 0;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--methods" && i + 1 < argc) methods = split(argv[++i], ',');
        else if (a == "--target" && i + 1 < argc) target = std::stod(argv[++i]);
        else if (a == "--h" && i + 1 < argc) h = std::stod(argv[++i]);
        else if (a == "--obj" && i + 1 < argc) objDir = argv[++i];
        else if (a == "--vtk" && i + 1 < argc) vtkDir = argv[++i];
        else if (a == "--limit" && i + 1 < argc) limit = std::stoi(argv[++i]);
        else if (a == "--no-layout") doLayout = false;
        else if (a == "--no-junction-domains") doJunctionDomains = false;
        else if (a == "--disk-templates") diskTemplates = true;
        else if (a == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--out DIR] [--methods SIPG,B1,B2] [--target EDGE] [--h SIZE]\n"
                         "       [--obj DIR] [--vtk DIR] [--limit N] [--no-layout]\n"
                         "       [--no-junction-domains] [--disk-templates] [model ...]\n";
            return 0;
        } else only.push_back(a);
    }
    ensureDir(outDir);
    if (!objDir.empty()) ensureDir(objDir);
    if (!vtkDir.empty()) ensureDir(vtkDir);
    Verdicts v;

    banner("E5  Multi-material domains");

    // The corpus, plus the five authored junction domains. The outline asks for
    // at least five domains before E5 gets a quantitative table; data/meshes/
    // multimat has fifteen, so the authored ones are here for Fig. 6's close-ups
    // and because a junction built to be one junction is the clearest place to
    // measure the junction residual.
    std::vector<Case> cases;
    for (const Model &mm : select(corpus("multimat"), only)) {
        std::string why;
        std::shared_ptr<Mesh> mesh = load(mm, why);
        if (!mesh) { v.warn(mm.name + ": " + why); continue; }
        if (!isMultiMaterial(*mesh)) {
            v.warn(mm.name + ": only one material id; skipped");
            continue;
        }
        cases.push_back({mm.name, mesh, "corpus"});
    }
    if (doJunctionDomains && only.empty()) {
        for (const Domain &d : junctionDomains(h)) {
            if (!isMultiMaterial(*d.mesh)) {
                v.warn(d.name + ": the mesher gave it one material; skipped");
                continue;
            }
            cases.push_back({d.name, d.mesh, d.note});
        }
    }
    if (limit > 0 && static_cast<int>(cases.size()) > limit) cases.resize(limit);
    if (cases.empty()) { std::cout << kFail << " no multi-material domains\n"; return 2; }

    std::cout << cases.size() << " domain(s), methods:";
    for (const std::string &m : methods) std::cout << " " << m;
    std::cout << "\n";

    // -----------------------------------------------------------------------
    // Parts 1 and 2 -- the field
    // -----------------------------------------------------------------------
    banner("E5  Part 1 and 2: interface alignment and the junction vertices");

    Csv fcsv(outDir + "/E5_field.csv",
             {"domain", "kind", "method", "vertices", "triangles", "materials",
              "interface_edges", "junctions", "junction_obliquity_max_deg",
              "interface_align_max_deg", "interface_align_p95_deg", "interface_align_mean_deg",
              "boundary_align_max_deg",
              "junction_residual_max_deg", "junction_residual_mean_deg",
              "junctions_over_1deg", "singularities", "ph_residual", "energy", "seconds"});

    struct FieldRow {
        std::string domain, method;
        double ifaceMax = 0, ifaceP95 = 0, junctionMax = 0, junctionMean = 0;
        int junctionsOver = 0, junctions = 0;
    };
    std::vector<FieldRow> frows;
    std::map<std::string, Eigen::VectorXcd> keptField;   // for the layout part

    for (const Case &c : cases) {
        const std::vector<int> ie = interfaceEdges(*c.mesh);
        const Junctions J = junctions(*c.mesh);
        std::set<int> mats(c.mesh->triangleMatId.begin(), c.mesh->triangleMatId.end());

        std::cout << "\n" << c.name << "  (" << c.mesh->vertices.size() << " v, "
                  << c.mesh->triangles.size() << " t, " << mats.size() << " materials, "
                  << ie.size() << " interface edges, " << J.vertices.size() << " junctions)"
                  << (c.note == "corpus" ? "" : "  -- " + c.note) << "\n";

        double obliq = 0.0;
        for (double x : J.obliquity) obliq = std::max(obliq, x);
        if (!J.vertices.empty())
            std::cout << "  worst junction obliquity " << num(obliq, 4)
                      << " deg (0 means every sector is a multiple of a right angle,"
                      << " so one cross satisfies them all)\n";

        Table t({"method", "iface align max", "p95", "mean", "dS align max",
                 "junction max", "junction mean", "junctions > 1 deg", "#sing", "PH", "energy"});

        for (const std::string &name : methods) {
            MethodOptions o;
            o.convergenceTol = kTol;
            const FieldRun run = runMethod(name, c.mesh, o);
            if (!run.ok) { v.warn(c.name + "/" + name + ": " + run.error); continue; }

            const metrics::AlignmentError ia = metrics::alignmentError(*c.mesh, run.u, ie);
            const metrics::AlignmentError ba = metrics::boundaryAlignmentError(*c.mesh, run.u);
            const metrics::SingularityReport sr = metrics::singularities(*c.mesh, run.u);
            const double energy = metrics::commonEnergy(
                *c.mesh, metrics::edgeWeights(*c.mesh, kGammaEval), run.u,
                interfaceEdgeFlags(*c.mesh));

            // The junction residual. B1 is measured on its *vertex* field, which
            // is the object the claim is about; a face method is measured on its
            // faces. Measuring B1's converted face field instead would hide the
            // very thing being shown, because the conversion has already
            // averaged the conflict away.
            std::vector<double> res;
            for (std::size_t k = 0; k < J.vertices.size(); ++k) {
                const double x = (run.hasConversion && run.p1Vertex.size() > 0)
                    ? junctionResidualVertex(*c.mesh, run.p1Vertex, J.vertices[k], J.edgesAt[k])
                    : junctionResidualFace(*c.mesh, run.u, J.edgesAt[k]);
                res.push_back(deg(x));
            }
            const double jMax = res.empty() ? 0.0 : *std::max_element(res.begin(), res.end());
            const double jMean = mean(res);
            int over = 0;
            for (double x : res) if (x > 1.0) ++over;

            t.row({name, num(deg(ia.maxAngle), 4), num(deg(ia.p95Angle), 4),
                   num(deg(ia.meanAngle), 4), num(deg(ba.maxAngle), 4),
                   num(jMax, 4), num(jMean, 4), num(over) + "/" + num((int)res.size()),
                   num((int)sr.interior.size()), num(sr.poincareHopfResidual),
                   num(energy, 6)});

            fcsv.row({{"domain", c.name}, {"kind", c.note == "corpus" ? "corpus" : "authored"},
                      {"method", name},
                      {"vertices", num((int)c.mesh->vertices.size())},
                      {"triangles", num((int)c.mesh->triangles.size())},
                      {"materials", num((int)mats.size())},
                      {"interface_edges", num((int)ie.size())},
                      {"junctions", num((int)J.vertices.size())},
                      {"junction_obliquity_max_deg", num(obliq, 6)},
                      {"interface_align_max_deg", num(deg(ia.maxAngle), 6)},
                      {"interface_align_p95_deg", num(deg(ia.p95Angle), 6)},
                      {"interface_align_mean_deg", num(deg(ia.meanAngle), 6)},
                      {"boundary_align_max_deg", num(deg(ba.maxAngle), 6)},
                      {"junction_residual_max_deg", num(jMax, 6)},
                      {"junction_residual_mean_deg", num(jMean, 6)},
                      {"junctions_over_1deg", num(over)},
                      {"singularities", num((int)sr.interior.size())},
                      {"ph_residual", num(sr.poincareHopfResidual)},
                      {"energy", num(energy, 10)},
                      {"seconds", num(run.totalSeconds, 4)}});

            frows.push_back({c.name, name, deg(ia.maxAngle), deg(ia.p95Angle),
                             jMax, jMean, over, (int)res.size()});

            if (name == methods.front()) keptField[c.name] = run.u;
            if (!vtkDir.empty())
                writeFieldVTK(vtkDir + "/E5_" + c.name + "_" + name + ".vtk", *c.mesh, run.u);

            if (name == "SIPG") {
                // The 95th percentile rather than the max, and the max reported
                // beside it. A triangle carrying two interface edges that are
                // not parallel -- a corner of the interface network, or a
                // junction -- is pinned to the kappa-weighted average of the
                // two tangents, so it cannot be exactly tangent to either and
                // the max picks that triangle out. It is a handful of edges on
                // a model and it is the same corner-averaging limitation
                // Sec. 7 lists for dS; the rest of the network is exact.
                int over = 0;
                for (int e : ie) {
                    const Point d = c.mesh->vertices[c.mesh->edges[e][1]]
                                  - c.mesh->vertices[c.mesh->edges[e][0]];
                    if (normP(d) < 1e-14) continue;
                    const double th = std::atan2(d[1], d[0]);
                    for (int side = 0; side < 2; ++side) {
                        const int t = c.mesh->edgeTriangles[e][side];
                        if (t >= 0 && deg(metrics::crossToDirection(run.u[t], th)) > 1.0) ++over;
                    }
                }
                v.check(deg(ia.p95Angle) < 1e-6,
                        c.name + "/SIPG: 95% of interface edges are a cross axis to rounding");
                if (over > 0)
                    v.warn(c.name + "/SIPG: " + num(over) + " of " + num(2 * (int)ie.size())
                           + " interface-edge sides are more than 1 degree off, worst "
                           + num(deg(ia.maxAngle), 3)
                           + " deg -- corner triangles carrying two non-parallel interface edges");
                v.check(sr.poincareHopfResidual == 0, c.name + "/SIPG: Poincare-Hopf holds");
            }
        }
        t.print();
    }

    heading("Part 2 summary: the junction residual");
    {
        Table t({"method", "domains", "junctions", "residual max (deg)", "residual mean (deg)",
                 "junctions over 1 deg"});
        for (const std::string &name : methods) {
            int domains = 0, junc = 0, over = 0;
            double worst = 0.0;
            std::vector<double> means;
            for (const FieldRow &r : frows) {
                if (r.method != name) continue;
                ++domains;
                junc += r.junctions;
                over += r.junctionsOver;
                worst = std::max(worst, r.junctionMax);
                means.push_back(r.junctionMean);
            }
            t.row({name, num(domains), num(junc), num(worst, 4), num(mean(means), 4),
                   num(over)});
        }
        t.print();
        std::cout << "\n  A face-based field answers the junction question per face and the\n"
                  << "  residual is rounding. A vertex-based field answers it once for every\n"
                  << "  incident material at the same vertex, and what is left is the\n"
                  << "  measurement Sec. 4.4 states as a proposition.\n"
                  << "\n  Read it with the obliquity column. A junction whose sectors are all\n"
                  << "  multiples of a right angle is a *control*: a cross is invariant under a\n"
                  << "  quarter turn, so one cross does satisfy every constraint there and a\n"
                  << "  vertex field has no difficulty with it. The proposition is about the\n"
                  << "  oblique junctions, and the residual should be read on those rows.\n"
                  << "\n  A second thing the tables above do not say for themselves: B1 has no\n"
                  << "  interface-alignment mechanism at all -- CrossField's only Dirichlet data\n"
                  << "  is dS. Its interface alignment error is therefore not a worse answer to\n"
                  << "  the same problem but the absence of the constraint, and its lower energy\n"
                  << "  on some domains is the energy of an unconstrained field. The comparison\n"
                  << "  that is fair is the junction residual, where both are asked what cross\n"
                  << "  they carry at a vertex the interfaces run through.\n";
    }

    // -----------------------------------------------------------------------
    // Part 3 -- the layout
    // -----------------------------------------------------------------------
    if (doLayout) {
        banner("E5  Part 3: the layout, and the mixed-zone count");

        Csv lcsv(outDir + "/E5_layout.csv",
                 {"domain", "method", "reached", "valid", "materials",
                  "patches", "mesh_quads", "mixed_quads", "pipeline_mixed_quads",
                  "irregular_interior", "irregular_boundary", "irregular_per_material",
                  "min_scaled_jacobian", "mean_scaled_jacobian", "inverted_quads",
                  "min_zone_dimension", "p01_zone_dimension", "seconds"});

        for (const Case &c : cases) {
            std::cout << "\n" << c.name << std::flush;
            Table t({"method", "reached", "valid", "patches", "quads", "mixed",
                     "irr int", "per material", "min SJ", "min zone", "s"});
            for (const std::string &name : methods) {
                MethodOptions fo;   // the shipped tolerance: this measures the pipeline
                const FieldRun run = runMethod(name, c.mesh, fo);
                if (!run.ok) { v.warn(c.name + "/" + name + ": " + run.error); continue; }

                LayoutOptions lo;
                lo.quadTargetEdge = target;
                lo.diskTemplates = diskTemplates;
                const LayoutResult L = runLayout(c.mesh, run.u, lo, !objDir.empty());
                const metrics::QuadMetrics &Q = L.quality;

                t.row({name, L.reachedStage, L.valid ? "yes" : "no", num(L.patches),
                       num(Q.quads), num(L.pipelineMixedQuads),
                       num(Q.irregularInterior), perMaterial(Q.irregularPerMaterial),
                       num(Q.minScaledJacobian, 4), num(Q.minZoneDimension, 4),
                       num(L.seconds, 2)});

                lcsv.row({{"domain", c.name}, {"method", name},
                          {"reached", L.reachedStage}, {"valid", num(L.valid)},
                          {"materials", num(L.pipelineMaterials)},
                          {"patches", num(L.patches)}, {"mesh_quads", num(Q.quads)},
                          {"mixed_quads", num(Q.mixedQuads)},
                          {"pipeline_mixed_quads", num(L.pipelineMixedQuads)},
                          {"irregular_interior", num(Q.irregularInterior)},
                          {"irregular_boundary", num(Q.irregularBoundary)},
                          {"irregular_per_material", "\"" + perMaterial(Q.irregularPerMaterial) + "\""},
                          {"min_scaled_jacobian", num(Q.minScaledJacobian, 6)},
                          {"mean_scaled_jacobian", num(Q.meanScaledJacobian, 6)},
                          {"inverted_quads", num(Q.invertedQuads)},
                          {"min_zone_dimension", num(Q.minZoneDimension, 6)},
                          {"p01_zone_dimension", num(Q.p01ZoneDimension, 6)},
                          {"seconds", num(L.seconds, 3)}});

                // R1 is asserted only where there *is* a layout to assert it
                // about. A model the pipeline did not finish -- bubbles without
                // Stage 0c's disk templates leaves more than half of itself
                // unmeshed -- has mixed elements because it has holes in it,
                // which is a different failure and is reported as one.
                if (name == "SIPG" && Q.quads > 0) {
                    if (L.valid && L.unmeshedPatches == 0)
                        v.check(L.pipelineMixedQuads == 0,
                                c.name + "/SIPG: no element straddles a material interface (R1)");
                    else if (L.pipelineMixedQuads > 0)
                        v.warn(c.name + "/SIPG: " + num(L.pipelineMixedQuads)
                               + " mixed element(s), on a layout that stopped at "
                               + L.reachedStage + " with " + num(L.unmeshedPatches)
                               + " unmeshed patch(es) -- R1 is not measurable here");
                }

                if (!objDir.empty() && !L.quadCells.empty())
                    writeQuadOBJ(objDir + "/E5_" + c.name + "_" + name + ".obj", L);
            }
            std::cout << "\n";
            t.print();
        }
    }

    // -----------------------------------------------------------------------
    // E5(a) -- the boundary-condition ablation
    // -----------------------------------------------------------------------
    banner("E5(a)  Hard-pinned against weak Nitsche boundary conditions");
    {
        Csv acsv(outDir + "/E5a_bc_ablation.csv",
                 {"domain", "boundary_conditions", "energy", "singularities", "index4_sum",
                  "ph_residual", "iterations", "converged",
                  "boundary_align_max_deg", "interface_align_max_deg", "seconds"});

        // Five domains, as the outline asks: the first five of whatever the run
        // was given, which on a full run is the first five of the corpus.
        std::size_t n = std::min<std::size_t>(5, cases.size());
        Table t({"domain", "BC", "energy", "#sing", "sum 4I", "PH", "iters",
                 "dS align max", "iface align max", "s"});
        for (std::size_t i = 0; i < n; ++i) {
            const Case &c = cases[i];
            const std::vector<int> ie = interfaceEdges(*c.mesh);
            for (int hard = 1; hard >= 0; --hard) {
                MethodOptions o;
                o.convergenceTol = kTol;
                o.hardBoundary = (hard != 0);
                const FieldRun run = runSIPG(c.mesh, o);
                if (!run.ok) { v.warn(c.name + " BC ablation: " + run.error); continue; }
                const metrics::SingularityReport sr = metrics::singularities(*c.mesh, run.u);
                const double energy = metrics::commonEnergy(
                    *c.mesh, metrics::edgeWeights(*c.mesh, kGammaEval), run.u,
                    interfaceEdgeFlags(*c.mesh));
                const metrics::AlignmentError ba = metrics::boundaryAlignmentError(*c.mesh, run.u);
                const metrics::AlignmentError ia = metrics::alignmentError(*c.mesh, run.u, ie);
                const std::string kind = hard ? "hard" : "weak";

                t.row({c.name, kind, num(energy, 6), num((int)sr.interior.size()),
                       num(sr.interiorIndex4Sum), num(sr.poincareHopfResidual),
                       num(run.iterations), num(deg(ba.maxAngle), 4),
                       num(deg(ia.maxAngle), 4), num(run.totalSeconds, 3)});
                acsv.row({{"domain", c.name}, {"boundary_conditions", kind},
                          {"energy", num(energy, 10)},
                          {"singularities", num((int)sr.interior.size())},
                          {"index4_sum", num(sr.interiorIndex4Sum)},
                          {"ph_residual", num(sr.poincareHopfResidual)},
                          {"iterations", num(run.iterations)},
                          {"converged", num(run.converged)},
                          {"boundary_align_max_deg", num(deg(ba.maxAngle), 6)},
                          {"interface_align_max_deg", num(deg(ia.maxAngle), 6)},
                          {"seconds", num(run.totalSeconds, 4)}});
            }
        }
        t.print();
        std::cout << "\n  The results in this paper use the hard pin. The weak form is the same\n"
                  << "  assembly with the row elimination left out, so the two differ only in\n"
                  << "  how far alignment is allowed to trade against smoothness at dS.\n";
    }

    banner("E5 summary");
    std::cout << (v.failures ? kFail : kPass) << " " << v.failures << " failure(s), "
              << v.warnings << " warning(s). CSVs in " << outDir << "/\n";
    return v.failures ? 1 : 0;
}
