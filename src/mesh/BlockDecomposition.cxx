#include "mesh/BlockDecomposition.hxx"

#include <fstream>

std::vector<Point> BlockDecomposition::sidePolyline(int block, int side) const {
    if (block < 0 || block >= static_cast<int>(blocks.size())) return {};
    if (side < 0 || side > 3) return {};
    const Block &b = blocks[block];
    const int e = b.edges[side];
    if (e < 0 || e >= static_cast<int>(edges.size())) return {};
    const MacroEdge &me = edges[e];
    if (!b.flip[side]) return me.points;
    return std::vector<Point>(me.points.rbegin(), me.points.rend());
}

bool BlockDecomposition::writeEdgesOBJ(const std::string &path) const {
    std::ofstream out(path);
    if (!out) return false;
    out.precision(17);
    out << "# " << source << " block decomposition: " << blocks.size() << " blocks, "
        << edges.size() << " macro edges, " << vertices.size() << " macrovertices\n";
    int base = 1;
    for (const MacroEdge &me : edges) {
        if (me.points.size() < 2) continue;
        for (const Point &p : me.points) out << "v " << p[0] << " " << p[1] << " 0\n";
        out << "l";
        for (size_t i = 0; i < me.points.size(); ++i) out << " " << (base + static_cast<int>(i));
        out << "\n";
        base += static_cast<int>(me.points.size());
    }
    return static_cast<bool>(out);
}

bool BlockDecomposition::writeBlocksOBJ(const std::string &path) const {
    std::ofstream out(path);
    if (!out) return false;
    out.precision(17);
    out << "# " << source << " blocks: " << blocks.size() << "\n";
    int base = 1;
    for (size_t b = 0; b < blocks.size(); ++b) {
        std::vector<Point> loop;
        for (int s = 0; s < 4; ++s) {
            const std::vector<Point> side = sidePolyline(static_cast<int>(b), s);
            for (size_t k = 0; k + 1 < side.size(); ++k) loop.push_back(side[k]);
        }
        if (loop.size() < 3) continue;
        for (const Point &p : loop) out << "v " << p[0] << " " << p[1] << " 0\n";
        out << "l";
        for (size_t i = 0; i < loop.size(); ++i) out << " " << (base + static_cast<int>(i));
        out << " " << base << "\n";
        base += static_cast<int>(loop.size());
    }
    return static_cast<bool>(out);
}
