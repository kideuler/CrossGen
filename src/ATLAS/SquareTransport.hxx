#ifndef __SQUARE_TRANSPORT_HXX__
#define __SQUARE_TRANSPORT_HXX__

#include <array>

// Integer square transport, Sec. 2.1 of
// docs/square_transport_2d_theory_and_implementation.md.
//
// Every carrier cell q has local coordinates u_q in [0,1]^2, with corner c of
// the cell at
//
//     (0,0), (1,0), (1,1), (0,1)   for c = 0, 1, 2, 3
//
// -- the same unit-square convention as mesh::QuadMesh and MFEM, so a cell's
// corner order *is* its coordinate system and nothing else has to be stored to
// define it. For two cells q, r sharing an edge, the transition
//
//     g_qr(u_r) = R_qr u_r + t_qr,    R_qr in C_4,  t_qr in Z^2
//
// maps r's coordinates into q's frame: it identifies the common edge and places
// r's abstract square on the far side of it. Both halves are integers, so a
// chain of them composes exactly and a development (Sec. 5.1) is a finite
// computation with a yes-or-no answer rather than a tolerance.
//
// ### Deriving g_qr from the two side indices
//
// Side k of a cell runs from corner k to corner k+1, and its direction in the
// cell's own frame is R^k e_1. If q's side i is glued to r's side j the two
// run in opposite directions (both cells are counter-clockwise), so R_qr must
// carry R^j e_1 onto -R^i e_1 = R^(i+2) e_1:
//
//     R_qr = R^(i + 2 - j)
//
// and the translation is whatever sends r's corner j onto q's corner i+1:
//
//     t_qr = corner_q(i+1) - R_qr corner_r(j).
//
// Then g_qr also sends r's corner j+1 onto q's corner i, and because a quarter
// turn preserves orientation r's square lands outside q's. transportAcross()
// is that formula and nothing more; the carrier stores its result per side so
// Stage 4 reads the table rather than re-deriving it (Sec. 12, "Transport
// table").
//
// ### Overflow
//
// Sec. 5.2 asks for checked integer arithmetic. A development's coordinates are
// bounded by the number of cells in it, which is bounded by the carrier -- a
// few hundred thousand at most -- so 64-bit arithmetic cannot overflow on any
// input this code will meet. RectangleCertifier still refuses a coordinate
// beyond 2^40 rather than trusting that argument silently.

typedef std::array<long long, 2> IPoint;

struct SquareTransport {
    int rot = 0;          // k in {0,1,2,3}: R^k with R = [[0,-1],[1,0]]
    long long tx = 0;
    long long ty = 0;

    SquareTransport() = default;
    SquareTransport(int k, long long x, long long y) : rot(k & 3), tx(x), ty(y) {}

    static IPoint rotate(const IPoint &p, int k) {
        switch (k & 3) {
            case 0: return p;
            case 1: return IPoint{-p[1], p[0]};
            case 2: return IPoint{-p[0], -p[1]};
            default: return IPoint{p[1], -p[0]};
        }
    }

    IPoint apply(const IPoint &p) const {
        const IPoint r = rotate(p, rot);
        return IPoint{r[0] + tx, r[1] + ty};
    }

    // (*this) o b: apply b first, then this. G_r = G_q o g_qr is written
    // G_q.compose(g_qr).
    SquareTransport compose(const SquareTransport &b) const {
        const IPoint t = rotate(IPoint{b.tx, b.ty}, rot);
        return SquareTransport(rot + b.rot, t[0] + tx, t[1] + ty);
    }

    // x = R^k y + t  <=>  y = R^-k x - R^-k t.
    SquareTransport inverse() const {
        const int k = (4 - rot) & 3;
        const IPoint t = rotate(IPoint{-tx, -ty}, k);
        return SquareTransport(k, t[0], t[1]);
    }

    static SquareTransport translation(long long x, long long y) {
        return SquareTransport(0, x, y);
    }

    bool isIdentity() const { return rot == 0 && tx == 0 && ty == 0; }
    bool operator==(const SquareTransport &o) const {
        return rot == o.rot && tx == o.tx && ty == o.ty;
    }
    bool operator!=(const SquareTransport &o) const { return !(*this == o); }
};

// Corner c of the unit square, in the cell's own integer coordinates.
inline IPoint squareCorner(int c) {
    switch (c & 3) {
        case 0: return IPoint{0, 0};
        case 1: return IPoint{1, 0};
        case 2: return IPoint{1, 1};
        default: return IPoint{0, 1};
    }
}

// g_qr for q's side `sideQ` glued to r's side `sideR`. See the header comment.
inline SquareTransport transportAcross(int sideQ, int sideR) {
    const int k = (sideQ + 2 - sideR) & 3;
    const IPoint cq = squareCorner(sideQ + 1);
    const IPoint cr = SquareTransport::rotate(squareCorner(sideR), k);
    return SquareTransport(k, cq[0] - cr[0], cq[1] - cr[1]);
}

#endif // __SQUARE_TRANSPORT_HXX__
