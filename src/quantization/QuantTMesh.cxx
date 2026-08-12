#include "quantization/QuantTMesh.hxx"

#include <Eigen/Sparse>

#include <sstream>

int QuantTMesh::addEdge(double xIdeal) {
    Edge e;
    e.xIdeal = xIdeal;
    edges.push_back(e);
    finalized_ = false;
    return static_cast<int>(edges.size()) - 1;
}

int QuantTMesh::addFace(std::vector<int> s0, std::vector<int> s1,
                        std::vector<int> s2, std::vector<int> s3) {
    Face f;
    f.sides = {std::move(s0), std::move(s1), std::move(s2), std::move(s3)};
    faces.push_back(std::move(f));
    finalized_ = false;
    return static_cast<int>(faces.size()) - 1;
}

bool QuantTMesh::finalize(std::string *error) {
    auto fail = [&](const std::string &msg) {
        if (error) *error = msg;
        finalized_ = false;
        return false;
    };

    rows.clear();
    for (Edge &e : edges) {
        e.row = {PHANTOM, PHANTOM};
        e.sign = {0, 0};
        if (!(e.xIdeal > 0.0)) return fail("edge with non-positive xIdeal");
    }

    for (size_t fi = 0; fi < faces.size(); ++fi) {
        const Face &f = faces[fi];
        std::vector<char> seen(edges.size(), 0);
        for (int s = 0; s < 4; ++s) {
            if (f.sides[s].empty()) {
                std::ostringstream os;
                os << "face " << fi << " has an empty side " << s;
                return fail(os.str());
            }
            for (int e : f.sides[s]) {
                if (e < 0 || e >= static_cast<int>(edges.size())) {
                    std::ostringstream os;
                    os << "face " << fi << " references unknown edge " << e;
                    return fail(os.str());
                }
                if (seen[e]) {
                    std::ostringstream os;
                    os << "edge " << e << " appears twice in face " << fi;
                    return fail(os.str());
                }
                seen[e] = 1;
            }
        }
        for (int axis = 0; axis < 2; ++axis) {
            Row r;
            r.pos = f.sides[axis];
            r.neg = f.sides[axis + 2];
            r.face = static_cast<int>(fi);
            r.axis = axis;
            rows.push_back(std::move(r));
        }
    }

    for (size_t ri = 0; ri < rows.size(); ++ri) {
        for (int pass = 0; pass < 2; ++pass) {
            for (int ei : (pass == 0 ? rows[ri].pos : rows[ri].neg)) {
                Edge &e = edges[ei];
                int slot = e.row[0] == PHANTOM ? 0 : e.row[1] == PHANTOM ? 1 : -1;
                if (slot < 0) {
                    std::ostringstream os;
                    os << "edge " << ei << " appears in more than two faces";
                    return fail(os.str());
                }
                e.row[slot] = static_cast<int>(ri);
                e.sign[slot] = pass == 0 ? +1 : -1;
            }
        }
    }

    finalized_ = true;
    if (error) error->clear();
    return true;
}

long long QuantTMesh::sideSum(int face, int side) const {
    long long sum = 0;
    for (int e : faces[face].sides[side]) sum += edges[e].x;
    return sum;
}

bool QuantTMesh::consistent() const {
    if (rows.empty()) return true;
    std::vector<Eigen::Triplet<double>> trip;
    for (size_t ri = 0; ri < rows.size(); ++ri) {
        for (int e : rows[ri].pos) trip.emplace_back(static_cast<int>(ri), e, 1.0);
        for (int e : rows[ri].neg) trip.emplace_back(static_cast<int>(ri), e, -1.0);
    }
    Eigen::SparseMatrix<double> A(static_cast<int>(rows.size()),
                                  static_cast<int>(edges.size()));
    A.setFromTriplets(trip.begin(), trip.end());
    Eigen::VectorXd x(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) x[i] = edges[i].x;
    return (A * x).lpNorm<Eigen::Infinity>() == 0.0;
}

double QuantTMesh::objective() const {
    double sum = 0.0;
    for (const Edge &e : edges) {
        const double d = e.x / e.xIdeal - 1.0;
        sum += d * d;
    }
    return sum;
}
