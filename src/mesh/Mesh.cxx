// Mesh.cxx
#include "Mesh.hxx"
#include "VertexTriangleCSR.hxx"


// Helper struct for hashing an undirected edge (min,max)
struct EdgeKey {
	int a;
	int b;
	bool operator==(const EdgeKey &other) const { return a == other.a && b == other.b; }
};

struct EdgeKeyHash {
	std::size_t operator()(const EdgeKey &k) const {
		return static_cast<std::size_t>(k.a) * 73856093u ^ static_cast<std::size_t>(k.b) * 19349663u;
	}
};

// Parse an OBJ index token like "i", "i/j", or "i/j/k" and return the vertex index (1-based in OBJ)
static bool parse_obj_index(const std::string &tok, int &vertexIndexOut) {
	if (tok.empty()) return false;
	// find first '/' if present
	std::size_t slashPos = tok.find('/');
	std::string viStr = (slashPos == std::string::npos) ? tok : tok.substr(0, slashPos);
	try {
		// OBJ indices can be negative (relative). We only support positive absolute indices here.
		int vi = std::stoi(viStr);
		vertexIndexOut = vi;
		return true;
	} catch (...) {
		return false;
	}
}

// Parse the material id out of an OBJ material name such as "mat3": the
// trailing run of digits. Names without one keep the id already in effect.
static bool parse_material_id(const std::string &name, int &matIdOut) {
	std::size_t end = name.size();
	while (end > 0 && std::isdigit(static_cast<unsigned char>(name[end - 1]))) --end;
	if (end == name.size()) return false;
	try {
		matIdOut = std::stoi(name.substr(end));
		return true;
	} catch (...) {
		return false;
	}
}

Mesh::Mesh(const std::vector<Point> &verts, const std::vector<Triangle> &tris)
	: Mesh(verts, tris, std::vector<int>{}) {}

Mesh::Mesh(const std::vector<Point> &verts, const std::vector<Triangle> &tris,
           const std::vector<int> &matIds)
	: vertices(verts), triangles(tris) {

	// A mesh with no material information given is a single-material mesh.
	if (matIds.size() == triangles.size()) triangleMatId = matIds;
	else triangleMatId.assign(triangles.size(), 1);
	
	// Prepare adjacency; initialize with -1 for boundaries
	triangleAdjacency.resize(triangles.size(), std::array<int,3>{-1, -1, -1});
	triangleEdges.resize(triangles.size(), std::array<int,3>{-1, -1, -1});

	// Map edges to the triangle and edge id, and also track edge index
	struct EdgeInfo { int tri; int edgeId; int edgeIdx; };
	std::unordered_map<EdgeKey, EdgeInfo, EdgeKeyHash> edgeMap;
	edgeMap.reserve(triangles.size() * 3);

	auto makeKey = [](int u, int v) -> EdgeKey {
		if (u < v) return EdgeKey{u, v};
		return EdgeKey{v, u};
	};

	for (int t = 0; t < static_cast<int>(triangles.size()); ++t) {
		const Triangle &tri = triangles[t];
		int v0 = tri[0], v1 = tri[1], v2 = tri[2];
		EdgeKey e01 = makeKey(v0, v1);
		EdgeKey e12 = makeKey(v1, v2);
		EdgeKey e20 = makeKey(v2, v0);

		// For each edge, check if seen; if seen, set adjacency both ways
		auto handleEdge = [&](const EdgeKey &ek, int localEdgeId, int va, int vb){
			auto it = edgeMap.find(ek);
			if (it == edgeMap.end()) {
				// New edge - add to edges list
				int edgeIdx = static_cast<int>(edges.size());
				edges.push_back(std::array<int,2>{ek.a, ek.b}); // stored as (min, max)
				edgeTriangles.push_back(std::array<int,2>{t, -1}); // first triangle, second TBD
				edgeMap.emplace(ek, EdgeInfo{t, localEdgeId, edgeIdx});
				triangleEdges[t][localEdgeId] = edgeIdx;
			} else {
				// Found neighboring triangle
				int ot = it->second.tri;
				int oedge = it->second.edgeId;
				int edgeIdx = it->second.edgeIdx;
				triangleAdjacency[t][localEdgeId] = ot;
				triangleAdjacency[ot][oedge] = t;
				triangleEdges[t][localEdgeId] = edgeIdx;
				edgeTriangles[edgeIdx][1] = t; // set second triangle
			}
		};

		handleEdge(e01, 0, v0, v1);
		handleEdge(e12, 1, v1, v2);
		handleEdge(e20, 2, v2, v0);
	}

	// Build boundary edge list and isBoundaryEdge vector
	isBoundaryEdge.assign(edges.size(), false);
	for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
		if (edgeTriangles[e][1] == -1) {
			isBoundaryEdge[e] = true;
			boundaryEdges.push_back(e);
		}
	}

	// Classify boundary edges per triangle. A triangle with exactly 1 boundary edge is a regular boundary triangle.
	// A triangle with 2 or more boundary edges is considered a corner triangle; store first two boundary edges.
	for (int t = 0; t < static_cast<int>(triangles.size()); ++t) {
		int bEdges[3];
		int nb = 0;
		for (int e = 0; e < 3; ++e) {
			if (triangleAdjacency[t][e] == -1) {
				bEdges[nb++] = e;
			}
		}
		if (nb == 1) {
			boundaryTriangles.push_back(std::array<int,2>{t, bEdges[0]});
		} else if (nb >= 2) {
			cornerTriangles.push_back(std::array<int,3>{t, bEdges[0], bEdges[1]});
		}
	}

	// Collect boundary vertices from boundary edges
	std::unordered_set<int> bVertsSet;
	for (int edgeIdx : boundaryEdges) {
		const auto &edge = edges[edgeIdx];
		bVertsSet.insert(edge[0]);
		bVertsSet.insert(edge[1]);
	}
	boundaryVertices.assign(bVertsSet.begin(), bVertsSet.end());
	std::sort(boundaryVertices.begin(), boundaryVertices.end());

	// Create boolean flag vector for boundary vertices
	isBoundaryVertex.assign(vertices.size(), false);
	for (int bv : boundaryVertices) {
		if (bv >= 0 && bv < static_cast<int>(isBoundaryVertex.size())) {
			isBoundaryVertex[bv] = true;
		}
	}
	
	// Build CSR mapping of vertex -> incident triangles (CCW order)
	vertexTriangles = VertexTriangleCSR::buildFromMesh(*this);
}

Mesh::Mesh(const std::string &filename) {
	std::ifstream in(filename);
	if (!in) {
		throw std::runtime_error("Failed to open OBJ file: " + filename);
	}

	std::string line;
	std::vector<Point> tempVertices;
	std::vector<Triangle> tempTriangles;
	std::vector<int> tempMatIds;
	int currentMatId = 1; // faces before any usemtl belong to material 1

	while (std::getline(in, line)) {
		// Trim leading spaces
		auto ltrim = [](std::string &s){ s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch){ return !std::isspace(ch); })); };
		ltrim(line);
		if (line.empty() || line[0] == '#') continue;

		std::istringstream iss(line);
		std::string tag;
		iss >> tag;
		if (tag == "v") {
			// vertex: v x y [z]
			double x = 0.0, y = 0.0;
			iss >> x >> y; // 2D mesh expects x,y; ignore optional z if present
			Point p{ x, y };
			tempVertices.push_back(p);
		} else if (tag == "usemtl") {
			std::string matName;
			iss >> matName;
			parse_material_id(matName, currentMatId);
		} else if (tag == "f") {
			// face: expect triangles. If more than 3 vertices, triangulate fan-wise
			std::vector<int> faceIndices;
			std::string tok;
			while (iss >> tok) {
				int vi;
				if (!parse_obj_index(tok, vi)) continue;
				faceIndices.push_back(vi);
			}

			// Support negative indices (relative to end)
			auto resolveIndex = [&](int idx) -> int {
				int n = static_cast<int>(tempVertices.size());
				if (idx > 0) return idx - 1; // OBJ is 1-based
				// negative index: -1 refers to last vertex
				return n + idx; // idx is negative
			};

			if (faceIndices.size() < 3) {
				continue; // ignore invalid faces
			}

			// Triangulate polygon using a fan: (0,i,i+1)
			for (size_t i = 1; i + 1 < faceIndices.size(); ++i) {
				int a = resolveIndex(faceIndices[0]);
				int b = resolveIndex(faceIndices[i]);
				int c = resolveIndex(faceIndices[i + 1]);
				const auto &pa = tempVertices[a];
				const auto &pb = tempVertices[b];
				const auto &pc = tempVertices[c];
				double A2 = (pb[0]-pa[0])*(pc[1]-pa[1]) - (pb[1]-pa[1])*(pc[0]-pa[0]);
				if (A2 < 0.0) std::swap(b, c);
				tempTriangles.push_back(Triangle{a,b,c});
				tempMatIds.push_back(currentMatId);
			}
		}
		// ignore other tags (vt, vn, etc.)
	}

	// Delegate to the main constructor via placement new
	// This is a common pattern to reuse constructor logic
	this->~Mesh();
	new (this) Mesh(tempVertices, tempTriangles, tempMatIds);
}

int Mesh::findTriangleContainingPoint(const Point &p) const {
	if (triangles.empty()) return -1;

	// Start from the middle triangle (often a good starting point for balanced meshes)
	int currentTri = static_cast<int>(triangles.size()) / 2;
	
	const int maxIter = static_cast<int>(triangles.size()) + 100; // prevent infinite loops
	
	for (int iter = 0; iter < maxIter; ++iter) {
		const Triangle &tri = triangles[currentTri];
		const Point &v0 = vertices[tri[0]];
		const Point &v1 = vertices[tri[1]];
		const Point &v2 = vertices[tri[2]];

		// Compute barycentric coordinates using signed areas
		Point v0v1 = v1 - v0;
		Point v0v2 = v2 - v0;
		Point v0p = p - v0;

		double denom = cross2(v0v1, v0v2);
		if (std::abs(denom) < 1e-30) {
			// Degenerate triangle, try a neighbor or move to next triangle
			for (int e = 0; e < 3; ++e) {
				if (triangleAdjacency[currentTri][e] >= 0) {
					currentTri = triangleAdjacency[currentTri][e];
					break;
				}
			}
			continue;
		}

		double l1 = cross2(v0p, v0v2) / denom; // weight for v1
		double l2 = cross2(v0v1, v0p) / denom; // weight for v2
		double l0 = 1.0 - l1 - l2;             // weight for v0

		// Tolerance for being "inside"
		const double eps = -1e-10;

		// Check if point is inside this triangle
		if (l0 >= eps && l1 >= eps && l2 >= eps) {
			return currentTri;
		}

		// Point is outside - walk toward it by crossing the edge with most negative barycentric coord
		// The edge opposite to vertex i is edge i (connecting vertices (i+1)%3 and (i+2)%3)
		int crossEdge = -1;
		double minBary = 0.0;

		if (l0 < minBary) {
			minBary = l0;
			crossEdge = 0; // edge opposite v0, between v1 and v2 (local edge 1)
		}
		if (l1 < minBary) {
			minBary = l1;
			crossEdge = 1; // edge opposite v1, between v2 and v0 (local edge 2)
		}
		if (l2 < minBary) {
			minBary = l2;
			crossEdge = 2; // edge opposite v2, between v0 and v1 (local edge 0)
		}

		// Map from "opposite vertex" to actual local edge index
		// Edge 0 connects v0-v1 (opposite v2)
		// Edge 1 connects v1-v2 (opposite v0)
		// Edge 2 connects v2-v0 (opposite v1)
		int localEdge;
		if (crossEdge == 0) localEdge = 1;      // opposite v0 -> edge v1-v2 -> local edge 1
		else if (crossEdge == 1) localEdge = 2; // opposite v1 -> edge v2-v0 -> local edge 2
		else localEdge = 0;                      // opposite v2 -> edge v0-v1 -> local edge 0

		int neighbor = triangleAdjacency[currentTri][localEdge];
		if (neighbor < 0) {
			// Hit boundary - the point may be outside the mesh or exactly on this triangle's edge
			// Do a final precise check with slightly larger tolerance
			const double boundaryEps = 1e-9;
			if (l0 >= -boundaryEps && l1 >= -boundaryEps && l2 >= -boundaryEps) {
				return currentTri;
			}
			// Point is outside the mesh
			return -1;
		}

		currentTri = neighbor;
	}

	// Fallback: exhaustive search if walking failed (should rarely happen)
	for (int t = 0; t < static_cast<int>(triangles.size()); ++t) {
		const Triangle &tri = triangles[t];
		const Point &v0 = vertices[tri[0]];
		const Point &v1 = vertices[tri[1]];
		const Point &v2 = vertices[tri[2]];

		Point v0v1 = v1 - v0;
		Point v0v2 = v2 - v0;
		Point v0p = p - v0;

		double denom = cross2(v0v1, v0v2);
		if (std::abs(denom) < 1e-30) continue;

		double l1 = cross2(v0p, v0v2) / denom;
		double l2 = cross2(v0v1, v0p) / denom;
		double l0 = 1.0 - l1 - l2;

		const double eps = -1e-9;
		if (l0 >= eps && l1 >= eps && l2 >= eps) {
			return t;
		}
	}

	return -1; // Point not found in any triangle
}

