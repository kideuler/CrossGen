#ifndef __REALISATION_HXX__
#define __REALISATION_HXX__

#include <string>
#include <vector>

#include "ATLAS/BlockCover.hxx"
#include "ATLAS/CoarseDomain.hxx"
#include "ATLAS/SquareCarrier.hxx"

// Carry a carrier found on a CoarseDomain back onto the input domain, as a new
// carrier over the input that SquareCarrier::validate() then certifies.
//
// ### Blocks, not coarse cells
//
// Given the coarse carrier's cover, what is carried over vertex for vertex is
// only its macro complex -- block corners and the coarse edges along block
// sides -- and each block is then one transfinite grid over its four fine
// sides. The fine carrier has the same cells either way (a block's grid is
// the union of its coarse cells' grids), but the geometry differs where it
// matters: a curve that bulges into the domain between two samples is spread
// across the whole block, where a thin coarse cell along it would have taken
// all of the bulge and folded. The coarse interior was only ever a proxy.
//
// This is Sec. 11.2 ("assigning mesh counts after blocking") and Sec. 11.3
// ("distinguish exact chart geometry from sampled quads") applied to the
// coarse carrier as a whole: every coarse cell becomes an n x m tensor grid,
// so the fine carrier has exactly the coarse one's irregular vertices, and its
// base complex -- and so the blocking Stages 4 and 5 extract from it -- is the
// coarse layout's.
//
// ### Counts
//
// Union-find over the coarse edges, joining the opposite sides of each coarse
// cell, gives the classes of Sec. 11.2 that must share one count. A class's
// count is the larger of what its length asks for at the target size and what
// its boundary edges need: every input vertex on dS must survive (Sec. 1.1's
// geometric preservation, which validate() audits), so a coarse boundary edge
// that spans k input segments needs at least k fine edges, and any more are
// collinear points inserted on the input's own segments -- the one change to
// dS Sec. 8.4 allows without a neighbour to reconcile.
//
// ### Geometry
//
// A coarse boundary vertex goes to its place on the input boundary
// (CoarseDomain::locate), snapped onto an input vertex when it is within
// Options::snap of a segment of one. Boundary sides follow the input polyline
// exactly; interior sides start straight. Each cell's interior is the
// discrete transfinite interpolation of its four sides (the construction
// CavityFill::grid uses), which is exact on the sides but can fold where a
// straight coarse edge meets a curve that bends into the cell. So the interior
// is then untangled and smoothed by TMOP (mesh::TMOP, shape metric 002) with
// every boundary node fixed, and the better of the two states by inverted
// cells, then worst scaled Jacobian, is kept. Nothing here certifies anything:
// the verdict is validate() on the carrier built from getQuads().
class Realisation {
public:
    struct Options {
        // Target fine edge length; <= 0 means the input's mean edge length.
        double size = 0.0;
        // Snap a boundary point onto an input vertex within this fraction of
        // the segment it lies on.
        double snap = 0.3;
        // TMOP sweeps after the interpolation; 0 leaves the interpolation.
        int smoothingSweeps = 100;
    };

    struct Report {
        bool built = false;
        std::string reason;
        int coarseCells = 0, classes = 0, blocks = 0;
        int vertices = 0, cells = 0;
        int boundarySplits = 0, interfaceSplits = 0, snapped = 0;
        int invertedInterpolated = 0, inverted = 0;
        int invertedHarmonic = -1;      // -1: the harmonic map was not needed
        bool harmonic = false;          // the harmonic state replaced the interpolation
        int localRepairs = 0;           // local harmonic re-solves that were kept
        double minScaledJacobianInterpolated = 0.0;
        double minScaledJacobian = 0.0, meanScaledJacobian = 0.0;
        bool smoothed = false;          // the TMOP state was kept
        int tmopSweeps = 0, untangleSweeps = 0;
        double seconds = 0.0;
    };

    // With a cover of `carrier`, each block is realised as one grid; without,
    // each coarse cell.
    Realisation(const CoarseDomain &coarse, const SquareCarrier &carrier, const BlockCover *cover,
                const Options &opts);

    bool built() const { return report_.built; }
    const Report &getReport() const { return report_; }
    const SquareCarrier::Quads &getQuads() const { return quads_; }
    // Coarse carrier vertex -> its vertex in getQuads().
    const std::vector<int> &vertexMap() const { return vertexMap_; }

private:
    SquareCarrier::Quads quads_;
    std::vector<int> vertexMap_;
    Report report_;
};

#endif // __REALISATION_HXX__
