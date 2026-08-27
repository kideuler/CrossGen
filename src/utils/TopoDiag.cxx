#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <memory>
#include <limits>

#include "MERIDIAN/MERIDIAN.hxx"
#include "MERIDIAN/Separatrices.hxx"
#include "MERIDIAN/SubdomainLabels.hxx"
#include "mesh/Mesh.hxx"

int main(int argc, char **argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s mesh.obj\n", argv[0]); return 1; }
    auto mesh = std::make_shared<Mesh>(std::string(argv[1]));

    MERIDIAN::Options o;
    MERIDIAN m(mesh, o);
    if (!m.run()) std::fprintf(stderr, "(pipeline reported failure, continuing)\n");

    const SubdomainLabels &lab = m.getLabels();
    const LayoutEnergy &lay = m.getLayout();
    const Separatrices &sep = m.getSeparatrices();
    const std::vector<Point> &uv = lay.getUV();

    // extent, same definition as both stages: diagonal of the image
    Point lo{1e300, 1e300}, hi{-1e300, -1e300};
    for (const Point &p : uv) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    const double extent = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    std::printf("extent %.6g   layout report maxTopoResidual %.6e\n\n",
                extent, lay.getReport().maxTopoResidual);

    const std::vector<int> &cv = m.getImmersion().getConeVertices();

    std::printf("Gamma_topo paths, residual re-measured on the FINAL layout:\n");
    std::printf("  #   fromCone  toCone  subs  seams   seedGap      |resid|/extent\n");
    for (size_t i = 0; i < lab.topoPaths().size(); ++i) {
        const auto &tp = lab.topoPaths()[i];
        const double r = std::fabs(lab.topoResidual(uv, tp)) / extent;
        std::printf("  %2zu  %8d  %6d  %4zu  %5d  %10.3e  %12.3e\n",
                    i, cv[tp.fromCone], cv[tp.toCone], tp.subs.size(),
                    tp.seamCrossings, tp.seedGap / extent, r);
    }

    std::printf("\nStage 7 curves (snap %.1e of extent):\n", sep.getReport().snapTolerance / extent);
    std::printf("  from    dir  steps  seams  end        gap/extent   nearestCone\n");
    static const char *kD[4] = {"+u", "+v", "-u", "-v"};
    for (const auto &c : sep.curves()) {
        const char *e = c.end == Separatrices::End::Cone ? "cone"
                      : c.end == Separatrices::End::Boundary ? "dS"
                      : c.end == Separatrices::End::Capped ? "cap"
                      : c.end == Separatrices::End::Cycle ? "cycle" : "stuck";
        const double g = (c.end == Separatrices::End::Cone) ? c.gap : c.nearestConeGap;
        std::printf("  %6d  %3s  %5zu  %5d  %-6s  %12.3e   %6d\n",
                    cv[c.cone], kD[c.dir], c.steps.size(), c.seamCrossings, e,
                    std::isfinite(g) ? g / extent : -1.0,
                    c.nearestCone >= 0 ? cv[c.nearestCone] : -1);
    }
    return 0;
}
