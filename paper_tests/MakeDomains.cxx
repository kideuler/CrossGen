// MakeDomains -- write the mechanism set to data/meshes/<set>/ as .obj files.
//
// The domains of `paper::mechanismDomains` are built in memory, like E1's
// canonical ones and E5's junction ones. Writing them out is what lets the
// corpus-driven experiments (E3, E4, and E5 with --set) run on them without
// each one growing its own copy of the constructors, and it is what makes the
// set inspectable: a reviewer can open the .obj.
//
// The set is deliberately *not* written into data/meshes/singlemat or
// data/meshes/multimat. Those are the corpus the paper's aggregate tables
// average over, and a corpus chosen after seeing the results is not evidence.
// These domains were chosen to test a stated mechanism, they sweep it, and each
// family contains the setting at which the prediction says the two methods
// agree; they belong in a set of their own and are reported as one.
//
// Usage:  MakeDomains [--out DIR] [--h SIZE] [--list]

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common/Corpus.hxx"
#include "common/Domains.hxx"
#include "common/Methods.hxx"
#include "common/Metrics.hxx"
#include "common/Report.hxx"

using namespace paper;

namespace {

// The .obj convention data/meshes uses and `Mesh(const std::string&)` reads:
// vertices, then faces grouped by material with a `usemtl matN` before each
// group. Faces before any usemtl are material 1.
bool writeMeshOBJ(const std::string &path, const Mesh &m, const std::string &comment) {
    std::ofstream f(path);
    if (!f) return false;
    f << "# " << comment << "\n";
    f.setf(std::ios::fixed);
    f.precision(17);
    for (const Point &p : m.vertices) f << "v " << p[0] << " " << p[1] << " 0\n";

    std::map<int, std::vector<int>> byMat;
    for (std::size_t t = 0; t < m.triangles.size(); ++t)
        byMat[m.triangleMatId.empty() ? 1 : m.triangleMatId[t]].push_back(static_cast<int>(t));
    for (const auto &[mat, tris] : byMat) {
        f << "usemtl mat" << mat << "\n";
        for (int t : tris)
            f << "f " << m.triangles[t][0] + 1 << " " << m.triangles[t][1] + 1
              << " " << m.triangles[t][2] + 1 << "\n";
    }
    return true;
}

// What the mechanism actually is, measured on the domain rather than asserted:
// how far each boundary corner's interior angle is from a multiple of 90
// degrees, and the same for the sectors at each vertex where interfaces meet.
// A vertex-based field must round both; a face-based one rounds neither.
struct Obliquity {
    int corners = 0;             // boundary vertices that turn at all
    double cornerWorst = 0.0;    // worst |angle - nearest multiple of 90|, deg
    double cornerMean = 0.0;
    int junctions = 0;
    double junctionWorst = 0.0;  // worst sector deviation at a junction, deg
};

Obliquity obliquity(const Mesh &m) {
    Obliquity o;
    const std::vector<double> ang = metrics::interiorAngles(m);
    double sum = 0.0;
    for (std::size_t v = 0; v < m.vertices.size(); ++v) {
        if (!m.isBoundaryVertex[v]) continue;
        const double a = ang[v] * 180.0 / M_PI;
        if (std::fabs(a - 180.0) < 1.0) continue;      // a straight run, not a corner
        const double q = a / 90.0;
        const double dev = std::fabs(q - std::round(q)) * 90.0;
        ++o.corners;
        sum += dev;
        o.cornerWorst = std::max(o.cornerWorst, dev);
    }
    if (o.corners) o.cornerMean = sum / o.corners;

    // Interface junctions: a vertex where three or more materials meet, or
    // where an interface lands on dS. The sectors are the rays of every
    // tangency constraint there.
    const std::vector<int> ie = interfaceEdges(m);
    std::map<int, std::vector<int>> at;
    for (int e : ie) { at[m.edges[e][0]].push_back(e); at[m.edges[e][1]].push_back(e); }
    for (const auto &[v, edges] : at) {
        std::set<int> mats;
        const auto &vt = m.vertexTriangles;
        for (int k = vt.rowPtr[v]; k < vt.rowPtr[v + 1]; ++k) mats.insert(m.triangleMatId[vt.colIdx[k]]);
        const bool onB = m.isBoundaryVertex[v];
        if (!(edges.size() >= 3 || mats.size() >= 3 || (onB && !edges.empty()))) continue;

        std::vector<int> rays = edges;
        if (onB) {
            const auto &vtx = m.vertexTriangles;
            for (int k = vtx.rowPtr[v]; k < vtx.rowPtr[v + 1]; ++k)
                for (int c = 0; c < 3; ++c) {
                    const int e = m.triangleEdges[vtx.colIdx[k]][c];
                    if (m.isBoundaryEdge[e] && (m.edges[e][0] == v || m.edges[e][1] == v))
                        rays.push_back(e);
                }
        }
        std::vector<double> th;
        for (int e : rays) {
            const int w = (m.edges[e][0] == v) ? m.edges[e][1] : m.edges[e][0];
            const Point d = m.vertices[w] - m.vertices[v];
            if (normP(d) > 1e-14) th.push_back(std::atan2(d[1], d[0]));
        }
        std::sort(th.begin(), th.end());
        th.erase(std::unique(th.begin(), th.end()), th.end());
        if (th.size() < 2) continue;
        double worst = 0.0;
        for (std::size_t i = 0; i < th.size(); ++i) {
            double sec = th[(i + 1) % th.size()] - th[i];
            if (sec <= 0.0) sec += 2.0 * M_PI;
            const double q = sec * 2.0 / M_PI;
            worst = std::max(worst, std::fabs(q - std::round(q)) * 90.0);
        }
        ++o.junctions;
        o.junctionWorst = std::max(o.junctionWorst, worst);
    }
    return o;
}

} // namespace

int main(int argc, char **argv) {
    std::string outDir = std::string(PAPER_MESH_DIR) + "/mechanism";
    double h = 0.03;
    bool listOnly = false;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--h" && i + 1 < argc) h = std::stod(argv[++i]);
        else if (a == "--list") listOnly = true;
        else if (a == "--help") {
            std::cout << "Usage: " << argv[0] << " [--out DIR] [--h SIZE] [--list]\n";
            return 0;
        }
    }

    banner("The mechanism set");
    std::cout <<
        "Domains built to test where the two discretisations must differ. The\n"
        "mechanism is the boundary data: a vertex field rounds each corner into one\n"
        "of four quarter-turn classes, which is exact at 90, 180 and 270 degrees and\n"
        "wrong by up to 22.5 degrees half way between; a face field pins each edge to\n"
        "its own tangent. `corner dev` and `junction dev` below are how far this\n"
        "domain puts that rounding from exact, measured on the mesh. A row with 0 in\n"
        "both is a control and the prediction there is that the methods agree.\n\n"
        "Mesh size h = " << h << ".\n";

    std::vector<Domain> doms;
    try {
        doms = mechanismDomains(h);
    } catch (const std::exception &e) {
        std::cout << kFail << " building the set: " << e.what() << "\n";
        return 2;
    }

    if (!listOnly) ensureDir(outDir);
    Verdicts v;
    Csv csv(listOnly ? "/dev/null" : outDir + "/_manifest.csv",
            {"domain", "vertices", "triangles", "materials", "chi", "interface_edges",
             "corners", "corner_dev_worst_deg", "corner_dev_mean_deg",
             "junctions", "junction_dev_worst_deg", "note"});

    Table t({"domain", "v", "t", "mat", "chi", "iface e", "corners", "corner dev",
             "junctions", "junction dev", "what it tests"});
    for (const Domain &d : doms) {
        const Mesh &m = *d.mesh;
        std::set<int> mats(m.triangleMatId.begin(), m.triangleMatId.end());
        if (mats.empty()) mats.insert(1);
        const Obliquity o = obliquity(m);
        const int chi = metrics::eulerCharacteristic(m);
        const std::size_t ie = interfaceEdges(m).size();

        t.row({d.name, num((int)m.vertices.size()), num((int)m.triangles.size()),
               num((int)mats.size()), num(chi), num((int)ie), num(o.corners),
               num(o.cornerWorst, 3), num(o.junctions), num(o.junctionWorst, 3), d.note});
        csv.row({{"domain", d.name}, {"vertices", num((int)m.vertices.size())},
                 {"triangles", num((int)m.triangles.size())},
                 {"materials", num((int)mats.size())}, {"chi", num(chi)},
                 {"interface_edges", num((int)ie)}, {"corners", num(o.corners)},
                 {"corner_dev_worst_deg", num(o.cornerWorst, 6)},
                 {"corner_dev_mean_deg", num(o.cornerMean, 6)},
                 {"junctions", num(o.junctions)},
                 {"junction_dev_worst_deg", num(o.junctionWorst, 6)},
                 {"note", "\"" + d.note + "\""}});

        if (!listOnly) {
            const std::string path = outDir + "/" + d.name + ".obj";
            v.check(writeMeshOBJ(path, m, "mechanism set, " + d.name + ": " + d.note),
                    "wrote " + path);
        }
    }
    t.print();

    if (!listOnly)
        std::cout << "\n  Run the experiments on them with `--set mechanism`.\n";
    std::cout << (v.failures ? kFail : kPass) << " " << v.failures << " failure(s).\n";
    return v.failures ? 1 : 0;
}
