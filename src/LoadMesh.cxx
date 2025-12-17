// Simple utility to load a mesh and print its public data
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "polyvec/PolyVectors.hxx"

static void printMesh(const Mesh &m) {
	std::cout << "Vertices (" << m.vertices.size() << ")\n";
	for (size_t i = 0; i < m.vertices.size(); ++i) {
		const auto &p = m.vertices[i];
		std::cout << "  v[" << i << "]: (" << std::setprecision(17) << p[0] << ", " << p[1] << ")\n";
	}
	std::cout << "Triangles (" << m.triangles.size() << ")\n";
	for (size_t i = 0; i < m.triangles.size(); ++i) {
		const auto &t = m.triangles[i];
		std::cout << "  f[" << i << "]: (" << t[0] << ", " << t[1] << ", " << t[2] << ")\n";
	}
	std::cout << "Adjacency (" << m.triangleAdjacency.size() << ")\n";
	for (size_t i = 0; i < m.triangleAdjacency.size(); ++i) {
		const auto &a = m.triangleAdjacency[i];
		std::cout << "  adj[" << i << "]: (" << a[0] << ", " << a[1] << ", " << a[2] << ")\n";
	}
	std::cout << "Boundary edges (" << m.boundaryTriangles.size() << ")\n";
	for (size_t i = 0; i < m.boundaryTriangles.size(); ++i) {
		const auto &b = m.boundaryTriangles[i];
		std::cout << "  boundary[" << i << "]: (tri=" << b[0] << ", edge=" << b[1] << ")\n";
	}
}

static void printPolyField(const PolyField &field, const Mesh &mesh) {
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        PolyCoeffs coeffs = field.getPolyCoeffsForTriangle(i);
        std::cout << "Triangle " << i << " Poly Coeffs: x_0 = " << coeffs[0] << ", x_1 = " << coeffs[1] << "\n";
    }
}

static void printFieldVectors(const PolyField &field, const Mesh &mesh) {
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        const auto &vec = field.field[i];
        std::cout << "Triangle " << i << " Field Vectors: u = (" << vec.u[0] << ", " << vec.u[1] << "), v = (" << vec.v[0] << ", " << vec.v[1] << ")\n";
    }
}

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <mesh.obj>\n";
		return 1;
	}
	std::string path = argv[1];
	try {
		Mesh m(path);

        PolyField field(m);

        field.solveForPolyCoeffs();
        printPolyField(field, m);

        field.convertToFieldVectors();
        printFieldVectors(field, m);
        field.writeVTK("output.vtk");
	} catch (const std::exception &e) {
		std::cerr << "Error: " << e.what() << "\n";
		return 2;
	}
	return 0;
}
