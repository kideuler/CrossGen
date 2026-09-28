// CrossGen: Mesh2D using Gmsh
// Main program: load .geo, normalize to [-1,1]^2, mesh with h=2/np, save to output dir.

#include <gmsh.h>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <filesystem>
#include <unordered_map>
#include <map>
#include <array>
#include <fstream>
#include <iomanip>

namespace fs = std::filesystem;

static void setUniformMeshSize(double h) {
	// Set uniform mesh size on all points in the current model
	std::vector<std::pair<int, int>> points;
	gmsh::model::getEntities(points, 0); // dim=0 for points
	if (points.empty()) {
		// If points are not yet available, synchronize geometry to create them
		try {
			gmsh::model::geo::synchronize();
		} catch (...) {
			// ignore
		}
		gmsh::model::getEntities(points, 0);
	}
	// Apply size directly to point entity tags (vectorpair expected)
	if (!points.empty()) gmsh::model::mesh::setSize(points, h);
}

// Physical Surface groups, captured as (physical tag -> member surface tags).
// gmsh::model::occ::synchronize() has a quirk in this Gmsh version: after an
// occ::translate/occ::dilate call (even an identity one), the following
// synchronize() silently drops every dim-2 physical group even though the
// surface tags themselves are unaffected. normalizeToUnitSquare() takes this
// snapshot before transforming so the caller can restore it afterward.
static std::vector<std::pair<int, std::vector<int>>> capturePhysicalSurfaceGroups() {
	std::vector<std::pair<int, std::vector<int>>> groups;
	std::vector<std::pair<int, int>> physDimTags;
	gmsh::model::getPhysicalGroups(physDimTags, 2);
	for (const auto &dimTag : physDimTags) {
		std::vector<int> entityTags;
		gmsh::model::getEntitiesForPhysicalGroup(dimTag.first, dimTag.second, entityTags);
		groups.emplace_back(dimTag.second, entityTags);
	}
	return groups;
}

// Reinstates the snapshot unconditionally rather than only when
// getPhysicalGroups() reports empty: Gmsh can leave a physical group's tag
// reserved internally (so re-adding it collides) even once it has stopped
// showing up in that query, so checking "is it gone" first is not reliable.
// Clearing whatever dim-2 groups are currently registered before re-adding
// the snapshot is the safe order either way.
static void restorePhysicalSurfaceGroups(const std::vector<std::pair<int, std::vector<int>>> &saved) {
	if (saved.empty()) return;
	std::vector<std::pair<int, int>> current;
	gmsh::model::getPhysicalGroups(current, 2);
	// The geo-kernel path (single-material and geo-kernel multimat files)
	// keeps its groups through normalizeToUnitSquare, so leave it alone: a
	// touch-and-restore cycle on physical groups the frontal-Delaunay 2D
	// mesher hasn't run yet perturbs its advancing-front order and produces
	// a different (still valid, but needlessly different) mesh. Only the
	// occ::synchronize() quirk this guards against drops the count.
	if (current.size() == saved.size()) return;

	// Once occ::synchronize() has dropped a group, getPhysicalGroups() no
	// longer lists its tag, but the tag stays reserved internally: adding a
	// physical group with that same tag still fails as "already exists". Ask
	// for the removal of exactly the (dim, tag) pairs the snapshot names,
	// regardless of what the (unreliable, post-drop) query reports.
	std::vector<std::pair<int, int>> toRemove;
	toRemove.reserve(saved.size());
	for (const auto &group : saved) toRemove.emplace_back(2, group.first);
	try {
		gmsh::model::removePhysicalGroups(toRemove);
	} catch (...) {
		// Nothing registered under these tags; fine, addPhysicalGroup below will succeed.
	}
	for (const auto &group : saved) {
		gmsh::model::addPhysicalGroup(2, group.second, group.first);
	}
}

static void normalizeToUnitSquare() {
	// Compute bounding box of the whole model and translate/scale to fit [-1,1]^2
	double xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
	gmsh::model::getBoundingBox(-1, -1, xmin, ymin, zmin, xmax, ymax, zmax);

	double cx = 0.5 * (xmin + xmax);
	double cy = 0.5 * (ymin + ymax);
	double width = xmax - xmin;
	double height = ymax - ymin;
	double maxDim = std::max(width, height);
	if (maxDim <= 0) return;

		// Gather all entities to transform
		std::vector<std::pair<int, int>> all;
		gmsh::model::getEntities(all);

		// Apply transforms for both geometry kernels when available
		double s = 2.0 / maxDim;
		try {
			gmsh::model::geo::translate(all, -cx, -cy, 0.0);
			gmsh::model::geo::dilate(all, 0.0, 0.0, 0.0, s, s, 1.0);
			gmsh::model::geo::synchronize();
		} catch (...) {
			// ignore if GEO kernel not in use
		}
		try {
			gmsh::model::occ::translate(all, -cx, -cy, 0.0);
			gmsh::model::occ::dilate(all, 0.0, 0.0, 0.0, s, s, 1.0);
			gmsh::model::occ::synchronize();
		} catch (...) {
			// ignore if OCC kernel not in use
		}
}

// Material id of a surface entity: the tag of the first physical group it
// belongs to, or 1 when the geometry declares no physical groups at all (the
// single-material case). Multimaterial .geo files select ids by writing
// `Physical Surface(<id>) = {...};` for each region.
static int materialIdOfSurface(int surfaceTag) {
	std::vector<int> physicalTags;
	try {
		gmsh::model::getPhysicalGroupsForEntity(2, surfaceTag, physicalTags);
	} catch (...) {
		return 1;
	}
	if (physicalTags.empty()) return 1;
	return physicalTags.front();
}

static void writeOBJ(const std::string &path) {
	// Every node gmsh made, in its order. That includes nodes no triangle
	// uses -- the centre point of every Circle() in a .geo is one -- and those
	// are left out below: mesh::Mesh keeps every `v` line as a vertex, and a
	// vertex in no triangle still changes what UMBER and ATLAS compute (py/
	// dataset.csv's MAMBO faces, 2026-09-27).
	std::vector<std::size_t> nodeTags;
	std::vector<double> nodeCoords;
	std::vector<double> nodeParam;
	gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParam);

	// Gather the triangles of every surface entity, keyed by material id, so
	// that each material comes out as one contiguous `usemtl mat<id>` block.
	// They hold node tags until it is known which nodes are written.
	std::map<int, std::vector<std::array<std::size_t, 3>>> facesByMaterial;

	std::vector<std::pair<int, int>> surfaces;
	gmsh::model::getEntities(surfaces, 2);
	for (const auto &surface : surfaces) {
		const int surfaceTag = surface.second;
		const int matId = materialIdOfSurface(surfaceTag);
		std::vector<std::array<std::size_t, 3>> &faces = facesByMaterial[matId];

		std::vector<int> types;
		std::vector<std::vector<std::size_t>> elementTags, elementNodeTags;
		gmsh::model::mesh::getElements(types, elementTags, elementNodeTags, 2, surfaceTag);

		// Handle triangles (3-node and 6-node) and degrade quads to triangles
		for (std::size_t k = 0; k < types.size(); ++k) {
			int et = types[k];
			const auto &nodes = elementNodeTags[k];
			if (et == 2) {
				// 3-node triangles
				for (std::size_t j = 0; j + 2 < nodes.size(); j += 3) {
					faces.push_back({nodes[j + 0], nodes[j + 1], nodes[j + 2]});
				}
			} else if (et == 9) {
				// 6-node triangles: use corner nodes (1,2,3)
				for (std::size_t j = 0; j + 5 < nodes.size(); j += 6) {
					faces.push_back({nodes[j + 0], nodes[j + 1], nodes[j + 2]});
				}
			} else if (et == 3) {
				// 4-node quads: split into two triangles (1,2,3) and (1,3,4)
				for (std::size_t j = 0; j + 3 < nodes.size(); j += 4) {
					faces.push_back({nodes[j + 0], nodes[j + 1], nodes[j + 2]});
					faces.push_back({nodes[j + 0], nodes[j + 2], nodes[j + 3]});
				}
			} else {
				// Other element types are ignored for OBJ output
			}
		}
	}

	// The nodes some triangle uses, numbered in gmsh's order, so that a mesh
	// without stray nodes comes out exactly as it always has.
	std::unordered_map<std::size_t, int> idx;   // node tag -> 1-based OBJ index
	for (const auto &entry : facesByMaterial)
		for (const auto &f : entry.second)
			for (std::size_t tag : f) idx[tag] = 0;
	std::vector<std::size_t> order;   // positions in nodeTags, in OBJ order
	for (std::size_t i = 0; i < nodeTags.size(); ++i) {
		auto it = idx.find(nodeTags[i]);
		if (it == idx.end() || it->second != 0) continue;
		order.push_back(i);
		it->second = static_cast<int>(order.size());
	}

	std::ofstream out(path);
	if (!out) {
		throw std::runtime_error("Failed to open OBJ file for writing: " + path);
	}
	out.setf(std::ios::fixed, std::ios::floatfield);
	out << std::setprecision(17);

	out << "# OBJ generated by CrossGen Mesh2Dgmsh\n";

	for (std::size_t i : order) {
		out << "v " << nodeCoords[3 * i + 0] << ' ' << nodeCoords[3 * i + 1] << ' '
		    << nodeCoords[3 * i + 2] << '\n';
	}

	// Write faces, one `usemtl` block per material id
	for (const auto &entry : facesByMaterial) {
		if (entry.second.empty()) continue;
		out << "usemtl mat" << entry.first << '\n';
		for (const auto &f : entry.second) {
			out << "f " << idx[f[0]] << ' ' << idx[f[1]] << ' ' << idx[f[2]] << '\n';
		}
	}

	out.close();
}

int main(int argc, char **argv) {
	if (argc < 4) {
		std::cerr << "Usage: Mesh2Dgmsh <input.geo> <output_dir> <np>\n";
		return 1;
	}

	const std::string inputGeo = argv[1];
	const std::string outputDir = argv[2];
	const int np = std::stoi(argv[3]);
	if (np <= 0) {
		std::cerr << "Error: np must be a positive integer.\n";
		return 1;
	}

	try {
		gmsh::initialize();
		gmsh::option::setNumber("General.Terminal", 1);

		// Load the .geo file as the current model
		gmsh::open(inputGeo);

		// Use GEO factory if applicable
		gmsh::model::geo::synchronize();

		// Merge geometry (if the .geo creates multiple components, open already loads them).
		// Apply normalization: translate and scale to fit within [-1,1]^2. This can
		// drop OCC-kernel Physical Surface groups (see normalizeToUnitSquare's
		// comment above capturePhysicalSurfaceGroups), so save and restore them.
		auto savedPhysicalGroups = capturePhysicalSurfaceGroups();
		normalizeToUnitSquare();
		restorePhysicalSurfaceGroups(savedPhysicalGroups);

		// Set characteristic length h = 2/np
		const double h = 2.0 / static_cast<double>(np);
		setUniformMeshSize(h);

		// Generate 2D mesh
		gmsh::model::mesh::generate(2);

		// Ensure output directory exists
		fs::create_directories(outputDir);

	// Compose output file path: same basename as input .geo, with .obj in outputDir
		fs::path inPath(inputGeo);
		std::string baseName = inPath.stem().string();
		fs::path outPath = fs::path(outputDir) / (baseName + ".obj");

	// Write OBJ manually from mesh data
	writeOBJ(outPath.string());

		gmsh::finalize();
	} catch (const std::exception &e) {
		std::cerr << "Gmsh error: " << e.what() << "\n";
		try { gmsh::finalize(); } catch (...) {}
		return 1;
	} catch (...) {
		std::cerr << "Unknown error during meshing.\n";
		try { gmsh::finalize(); } catch (...) {}
		return 1;
	}

	return 0;
}

