#ifndef _MEDIAL_AXIS_MAP_HXX_
#define _MEDIAL_AXIS_MAP_HXX_

#include "medialaxis/MedialAxis.hxx"

#include <memory>
#include <vector>

// ─── Monotone boundary -> medial axis map ───────────────────────────────────
//
// The map phi : dD -> M(D) of "Planar quad meshing with guarantees", Sec. 3.
//
// The Voronoi cell of a boundary vertex v covers the boundary sub-polyline
// mid(prev, v) -> v -> mid(v, next) and, inside the domain, a chain of medial
// vertices: the circumcenters of the dual cells incident to v, ordered from
// the boundary edge (prev, v) round to (v, next). Exactly two Voronoi edges of
// the cell cross the boundary — the ones dual to those two boundary edges —
// and their intersection point c (possibly at infinity) is the center of the
// polar parameterization that carries the sub-polyline onto the chain
// (Sec. 3.3, Fig. 7). Both convex and concave vertices are covered: from c the
// chain and the sub-polyline both appear convex, so the map is a bijection.
//
// The map collapses only where the chain degenerates to a single medial
// vertex, i.e. where v's fan is a single dual cell. That is the discrete
// counterpart of a continuous medial axis endpoint: the whole sub-polyline
// maps to the one vertex, forming a polar section of the map (Fig. 14).
// Together these give monotonicity in the sense of Sec. 3.2 — restricted to
// each sub-polyline, phi is a bijection or has a single-point image.

// The restriction of phi to the Voronoi cell of one boundary vertex.
struct MedialFan {
    int boundaryVertex = -1;

    // Active dual cells incident to the boundary vertex, ordered from the
    // boundary edge (prev, v) to (v, next); `chain` holds their medial
    // vertices in the same order. Consecutive fans share their end vertices:
    // chain.back() of one fan is chain.front() of the next.
    std::vector<int> cells;
    std::vector<int> chain;

    // Where the two boundary-crossing Voronoi edges cross: the midpoints of
    // the two boundary edges at v. They delimit the sub-polyline this fan
    // maps, and they are the pre-images of the chain's end vertices.
    Point midPrev{0.0, 0.0};
    Point midNext{0.0, 0.0};

    // Polar center c. When the two Voronoi rays are (nearly) parallel the
    // center is at infinity and `parallelDir` carries the common ray
    // direction instead; the polar pencil degenerates to parallel lines.
    Point  center{0.0, 0.0};
    bool   centerFinite = false;
    Point  parallelDir{0.0, 0.0};

    // phi(v), the image of the boundary vertex itself.
    Point image{0.0, 0.0};

    // Pre-image on the sub-polyline of every chain vertex, same size as
    // `chain`; front() == midPrev and back() == midNext by construction. For
    // a collapsed fan it holds the single entry v.
    std::vector<Point> footpoints;

    // True when the whole sub-polyline maps to a single medial vertex.
    bool collapsed = false;
};

// One drawable medial radius: a boundary point joined to its image under phi.
// Emitted so that every radius appears exactly once: each fan contributes its
// midPrev radius, its interior footpoints and its vertex radius, and leaves
// the midNext radius to the following fan.
struct MedialSpoke {
    Point boundary{0.0, 0.0};
    Point medial{0.0, 0.0};
    bool collapsed = false;   // part of a polar section of the map
};

struct MedialAxisMapStats {
    int fans = 0;                // boundary vertices carrying a fan
    int collapsedFans = 0;       // fans whose chain is a single medial vertex
    int skippedVertices = 0;     // boundary vertices with unusable links/cells
    int fallbackProjections = 0; // polar rays that missed and fell back to a
                                 // closest-point projection
};

class MedialAxisMap {
public:
    // Builds phi for the axis as it currently stands, so run any
    // simplification (deduplicateMedialVertices) before constructing the map.
    explicit MedialAxisMap(std::shared_ptr<MedialAxis> axis);

    // phi evaluated at p, which should lie on the sub-polyline of
    // `boundaryVertex`. Returns p unchanged when that vertex has no fan.
    Point mapToAxis(int boundaryVertex, const Point &p) const;

    const std::vector<MedialFan> &fans() const { return fans_; }

    // Fan of a boundary vertex, or nullptr where none was built.
    const MedialFan *fanOf(int boundaryVertex) const;

    const std::vector<MedialSpoke> &spokes() const { return spokes_; }

    const MedialAxisMapStats &stats() const { return stats_; }

    std::shared_ptr<MedialAxis> axis;

private:
    void buildFan(int v);
    void buildSpokes();

    // The polar map on one fan. Sets `missed` instead of failing when no
    // chain segment lies on the polar line through p, in which case the
    // closest chain point stands in.
    Point mapOnFan(const MedialFan &fan, const Point &p, bool &missed) const;

    // Inverse direction: the boundary point of the fan's sub-polyline on the
    // polar line through medial point m.
    Point footpointOnBoundary(const MedialFan &fan, const Point &m,
                              const Point &vertexPos, bool &missed) const;

    std::vector<int> fanOfVertex_;   // fan index per mesh vertex, -1 if none
    std::vector<MedialFan> fans_;
    std::vector<MedialSpoke> spokes_;
    MedialAxisMapStats stats_;
};

#endif // _MEDIAL_AXIS_MAP_HXX_
