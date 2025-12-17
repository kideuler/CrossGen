// Mesh.cxx
#include "Mesh.hxx"


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

Mesh::Mesh(const std::string &filename) {
	// Reset containers
	vertices.clear();
	triangles.clear();
	triangleAdjacency.clear();
	boundaryTriangles.clear();
	cornerTriangles.clear();

	std::ifstream in(filename);
	if (!in) {
		throw std::runtime_error("Failed to open OBJ file: " + filename);
	}

	std::string line;
	std::vector<Point> tempVertices; // accumulate to allow negative indices if needed later

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
				triangles.push_back(Triangle{ a, b, c });
			}
		}
		// ignore other tags (vt, vn, etc.)
	}

	// Move vertices from temp to the public container
	vertices = std::move(tempVertices);

	// Prepare adjacency; initialize with -1 for boundaries
	triangleAdjacency.resize(triangles.size(), std::array<int,3>{-1, -1, -1});

	// Map edges to the triangle and edge id
	struct EdgeInfo { int tri; int edgeId; };
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
		auto handleEdge = [&](const EdgeKey &ek, int edgeId){
			auto it = edgeMap.find(ek);
			if (it == edgeMap.end()) {
				edgeMap.emplace(ek, EdgeInfo{t, edgeId});
			} else {
				// Found neighboring triangle
				int ot = it->second.tri;
				int oedge = it->second.edgeId;
				triangleAdjacency[t][edgeId] = ot;
				triangleAdjacency[ot][oedge] = t;
			}
		};

		handleEdge(e01, 0);
		handleEdge(e12, 1);
		handleEdge(e20, 2);
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
}

