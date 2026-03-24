#ifndef _MEDIAL_AXIS_HXX_
#define _MEDIAL_AXIS_HXX_

#include "mesh/Mesh.hxx"

class MedialAxis {
    public:
        std::shared_ptr<Mesh> mesh;
        std::vector<Point> medialVertices; // List of medial axis vertices
        std::vector<Edge> medialEdges; // List of medial axis edges

        MedialAxis(std::shared_ptr<Mesh> mesh);
    private:
        Point computeCircumcenter(int triIndex);
        
};

#endif // _MEDIAL_AXIS_HXX_