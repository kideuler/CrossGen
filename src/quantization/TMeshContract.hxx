#ifndef _TMESH_CONTRACT_HXX_
#define _TMESH_CONTRACT_HXX_

#include "quantization/QuantTMeshConvert.hxx"

// ─── Removing zero-length edges from a quantized T-mesh ─────────────────────
//
// Some T-meshes admit no strictly positive quantization: the consistency
// system itself forces certain edges to zero, whatever the algorithm. A
// cell with a wholly zero side is then not a cell at all -- it has no
// quads -- and leaves both a hole in the picture and a degenerate entry in
// the structure.
//
// This is the cleanup QGP recommends for exactly that situation (Sec.
// 7.1.1, citing the partition simplification of Myles et al. 2014): remove
// the zero-cells from the T-mesh instead of carrying them. A collapsed
// cell has zero width in one direction, so the two sides that survive it
// are the same line in parameter space; identifying them deletes the cell
// and makes its two neighbours adjacent. Their shared edge is given the
// mid-curve between them, so between them they cover the ground the
// removed cell stood on and no hole is left.
//
// The two surviving sides carry the same number of quads but need not be
// cut into edges the same way, so an edge is occasionally split first to
// give both the same breaks.
//
// What comes out is a T-mesh whose edges are all >= 1, still consistent,
// and still quantized by the values already computed: the surviving
// entries of the old solution are a strictly positive solution of the new
// system, so nothing needs re-solving.

struct ContractReport {
    int mergedCells = 0;    // collapsed cells removed by identifying their sides
    int pointCells = 0;     // cells that collapsed to a point in both directions
    int splitEdges = 0;     // edges cut so two sides could be identified
    int removedEdges = 0;   // zero-length edges gone from the structure
    int remainingZero = 0;  // zero edges the cleanup could not remove
    bool ok = false;
    std::string error;
};

// Rewrites `bq` in place. Safe to call on a T-mesh that has no zero edges,
// where it does nothing. The tmesh must be quantized and consistent.
ContractReport contractZeroEdges(BlockQuant &bq);

#endif  // _TMESH_CONTRACT_HXX_
