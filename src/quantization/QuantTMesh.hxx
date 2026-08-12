#ifndef _QUANT_TMESH_HXX_
#define _QUANT_TMESH_HXX_

#include <array>
#include <string>
#include <vector>

// ─── Generic T-mesh for quantization ────────────────────────────────────────
//
// The combinatorial input of Campen, Bommes & Kobbelt, "Quantized Global
// Parametrization" (SIGGRAPH Asia 2015), Sec. 5, stripped of everything the
// planar setting makes unnecessary: no transitions, no matchings, no
// geometric embedding. A T-mesh here is a set of edges and a set of
// four-sided faces, each face listing the edges of its four sides in cyclic
// order, so sides[0] is opposite sides[2] and sides[1] opposite sides[3]. A
// side carrying T-junctions is simply several edges.
//
// Each edge carries one differential parameter x -- its integer parametric
// length, the unknown of the quantization -- and its ideal length xIdeal,
// read off whatever continuous geometry the T-mesh came from (or 1 for a
// pure block decomposition).
//
// finalize() derives the consistency conditions of Sec. 5.3: per face, the
// lengths of opposite sides sum to the same value. That is two rows of the
// homogeneous linear Diophantine system Ax = 0, x >= 0, with coefficients
// in {-1, 0, +1} and every edge appearing in exactly two rows -- or one,
// when the edge lies on the domain boundary, in which case its second row
// slot stays PHANTOM: the outer face never needs balancing. The quantizer
// rests entirely on this two-rows-per-edge structure.
//
// The structure is deliberately minimal so that richer layouts convert down
// to it easily; see QuantTMeshConvert.hxx for QuadLayout and MedialAxisTMesh.

class QuantTMesh {
public:
    // Row index meaning "no row": the edge borders the outer face there.
    static constexpr int PHANTOM = -1;

    struct Edge {
        double xIdeal = 1.0;
        int x = 0;  // the quantized length; the result of TMeshQuantizer

        // The (at most two) rows this edge appears in, with its coefficient
        // there. Filled by finalize(); row[1] is PHANTOM for boundary edges,
        // and both are PHANTOM for an edge no face references.
        std::array<int, 2> row{PHANTOM, PHANTOM};
        std::array<int, 2> sign{0, 0};

        bool onBoundary() const { return row[1] == PHANTOM; }
    };

    struct Face {
        // Edge ids, ordered along each side; sides[i] and sides[i+2] are the
        // opposite pairs.
        std::array<std::vector<int>, 4> sides;
    };

    // One consistency condition: sum of x over `pos` == sum over `neg`.
    struct Row {
        std::vector<int> pos, neg;
        int face = -1;
        int axis = 0;  // 0 balances sides 0 and 2, 1 balances sides 1 and 3
    };

    int addEdge(double xIdeal = 1.0);
    int addFace(std::vector<int> s0, std::vector<int> s1,
                std::vector<int> s2, std::vector<int> s3);

    // Build the rows and the per-edge row references, and validate the
    // structure: four non-empty sides, known edge ids, no edge twice in one
    // face, no edge in more than two faces, positive ideals. Returns false
    // and describes the offence in `error` if the input is not a T-mesh.
    bool finalize(std::string *error = nullptr);
    bool finalized() const { return finalized_; }

    long long sideSum(int face, int side) const;

    // Whether the current x satisfies Ax = 0 (checked with Eigen, so it
    // really is the matrix statement and not a re-derivation of the rows).
    bool consistent() const;

    // Sum over edges of (x / xIdeal - 1)^2, the Sec. 6.3 quality measure.
    double objective() const;

    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<Row> rows;

private:
    bool finalized_ = false;
};

#endif  // _QUANT_TMESH_HXX_
