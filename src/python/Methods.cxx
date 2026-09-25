#include "python/Methods.hxx"

#include <cmath>

namespace pycg {

PyObject *unknownStage(const char *method, const std::string &stage) {
    PyErr_Format(PyExc_ValueError,
                 "crossgen.options(): '%s' has no stage '%s' (expected 'method', 'mesh' or "
                 "'smooth')",
                 method, stage.c_str());
    return nullptr;
}

double decompositionCoverage(const BlockDecomposition &D, const Mesh &mesh) {
    double modelArea = 0.0;
    for (const Triangle &t : mesh.triangles) {
        const Point &a = mesh.vertices[t[0]], &b = mesh.vertices[t[1]], &c = mesh.vertices[t[2]];
        modelArea += 0.5 * std::fabs(cross2(b - a, c - a));
    }
    if (!(modelArea > 0.0)) return 0.0;

    // The shoelace over the outline, side after side. Each side repeats the
    // corner the one before it ended on, and a repeated point adds nothing to
    // the sum, so the four polylines can be walked as they are.
    double blockArea = 0.0;
    for (int b = 0; b < static_cast<int>(D.blocks.size()); ++b) {
        double twice = 0.0;
        for (int s = 0; s < 4; ++s) {
            const std::vector<Point> side = D.sidePolyline(b, s);
            for (size_t i = 0; i + 1 < side.size(); ++i) twice += cross2(side[i], side[i + 1]);
        }
        blockArea += 0.5 * std::fabs(twice);
    }
    return blockArea / modelArea;
}

}  // namespace pycg
