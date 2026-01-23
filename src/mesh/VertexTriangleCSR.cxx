#include "VertexTriangleCSR.hxx"
#include "Mesh.hxx"

static inline double angleAround(const Point &center, const Point &p) {
    double dx = p[0] - center[0];
    double dy = p[1] - center[1];
    return std::atan2(dy, dx);
}

static inline Point triCentroid(const Mesh &m, const Triangle &t) {
    const Point &a = m.vertices[t[0]];
    const Point &b = m.vertices[t[1]];
    const Point &c = m.vertices[t[2]];
    return Point{ (a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0 };
}

VertexTriangleCSR VertexTriangleCSR::buildFromMesh(const Mesh &mesh) {
    VertexTriangleCSR csr;
    int nV = static_cast<int>(mesh.vertices.size());
    int nT = static_cast<int>(mesh.triangles.size());
    csr.rowPtr.resize(nV + 1, 0);

    // Gather incident triangles per vertex
    std::vector<std::vector<int>> incident(nV);
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        incident[tri[0]].push_back(t);
        incident[tri[1]].push_back(t);
        incident[tri[2]].push_back(t);
    }

    // Compute rowPtr sizes first
    for (int v = 0; v < nV; ++v) {
        csr.rowPtr[v+1] = csr.rowPtr[v] + static_cast<int>(incident[v].size());
    }
    csr.colIdx.resize(csr.rowPtr[nV]);

    // For each vertex, sort incident triangles CCW by triangle centroid angle
    for (int v = 0; v < nV; ++v) {
        const Point &pv = mesh.vertices[v];
        auto &vec = incident[v];
        std::vector<std::pair<double,int>> angTri;
        angTri.reserve(vec.size());
        for (int t : vec) {
            Point c = triCentroid(mesh, mesh.triangles[t]);
            double a = angleAround(pv, c);
            angTri.emplace_back(a, t);
        }
        std::sort(angTri.begin(), angTri.end(), [](const auto &A, const auto &B){ return A.first < B.first; });
        int start = csr.rowPtr[v];
        for (size_t i = 0; i < angTri.size(); ++i) {
            csr.colIdx[start + static_cast<int>(i)] = angTri[i].second;
        }
    }

    return csr;
}
