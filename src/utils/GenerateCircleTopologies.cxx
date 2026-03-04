#include "tracing/SeparatrixTrace.hxx"
#include "TestHelper.hxx"

// for file writing
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>

// Write the mesh as a VTK Unstructured Grid (.vtu) file
void writeMeshVTU(const std::string &filename, const Mesh &mesh) {
    std::ofstream out(filename);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << mesh.vertices.size()
        << "\" NumberOfCells=\"" << mesh.triangles.size() << "\">\n";

    // Points
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &v : mesh.vertices) {
        out << "          " << v[0] << " " << v[1] << " 0.0\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // Cells (triangles)
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto &tri : mesh.triangles) {
        out << "          " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 1; i <= mesh.triangles.size(); ++i) {
        out << "          " << i * 3 << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        out << "          5\n"; // VTK_TRIANGLE = 5
    }
    out << "        </DataArray>\n";
    out << "      </Cells>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

// Write singularity points as a VTK PolyData (.vtp) file
void writeSingularitiesVTP(const std::string &filename,
                           const std::vector<Singularity> &singularities) {
    std::ofstream out(filename);
    int n = static_cast<int>(singularities.size());

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <PolyData>\n";
    out << "    <Piece NumberOfPoints=\"" << n
        << "\" NumberOfVerts=\"" << n << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    // Points
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &s : singularities) {
        out << "          " << s.coordinates[0] << " " << s.coordinates[1] << " 0.0\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // Vertex cells (one per point so they are visible)
    out << "      <Verts>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < n; ++i) {
        out << "          " << i << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 1; i <= n; ++i) {
        out << "          " << i << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Verts>\n";

    // Point data: singularity index
    out << "      <PointData Scalars=\"index\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"index\" format=\"ascii\">\n";
    for (const auto &s : singularities) {
        out << "          " << s.singularityIndex << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </PointData>\n";

    out << "    </Piece>\n";
    out << "  </PolyData>\n";
    out << "</VTKFile>\n";
}

// Write a VTK multiblock (.vtm) file referencing a mesh and singularities file
void writeVTM(const std::string &filename,
              const std::string &meshFile,
              const std::string &singFile) {
    std::ofstream out(filename);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n";
    out << "  <vtkMultiBlockDataSet>\n";
    out << "    <DataSet index=\"0\" name=\"mesh\" file=\"" << meshFile << "\"/>\n";
    out << "    <DataSet index=\"1\" name=\"singularities\" file=\"" << singFile << "\"/>\n";
    out << "  </vtkMultiBlockDataSet>\n";
    out << "</VTKFile>\n";
}

int main(int argc, char **argv) {
    // take 4 arguments: start_angle, end_angle, angle_step, meshsize
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <start_angle (in degrees)> <end_angle (in degrees)> <angle_step (in degrees)> <meshsize>\n";
        return 1;
    }

    double start_angle = std::stod(argv[1]);
    double end_angle = std::stod(argv[2]);
    double angle_step = std::stod(argv[3]);
    double meshsize = std::stod(argv[4]);

    // check inputs
    if (start_angle < 0 || start_angle >= 360.0) {
        std::cerr << "Error: start_angle must be in [0, 360)\n";
        return 1;
    }
    if (end_angle <= start_angle || end_angle > 360.0) {
        std::cerr << "Error: end_angle must be in (start_angle, 360]\n";
        return 1;
    }
    if (angle_step <= 0 || angle_step >= (end_angle - start_angle)) {
        std::cerr << "Error: angle_step must be in (0, end_angle - start_angle)\n";
        return 1;
    }
    if (meshsize <= 0) {
        std::cerr << "Error: meshsize must be positive\n";
        return 1;
    }

    // open file to write results
    std::ofstream out("circle_topologies.csv");
    if (!out.is_open()) {
        std::cerr << "Error: Could not open output file for writing\n";
        return 1;
    }

    // write CSV header
    out << "angle,theta1,r1,theta2,r2,theta3,r3,theta4,r4\n";

    // create output directory for VTK files
    std::string vtkDir = "circle_topologies_vtk";
    mkdir(vtkDir.c_str(), 0755);

    // collect PVD entries (angle, vtm filename)
    std::vector<std::pair<double, std::string>> pvdEntries;

    // Create a mesh using the provided parameters
    int frameIdx = 0;
    for (double angle = start_angle; angle <= end_angle; angle += angle_step, ++frameIdx) {
        std::cout << "Generating mesh for angle: " << angle << " degrees\n";

        // generate mesh
        auto mesh = TestHelper::createEllipse(
            0.0, 0.0,  // Center of the ellipse
            1.0, 1.0,  // Semi-major and semi-minor axes
            angle,    // Full sweep
            meshsize   // Mesh size
        );

        // create crossfield with mbo method
        auto crossField = std::make_shared<CrossField>(mesh, 1000);
        crossField->initialize(1); // method 1 = MBO
        crossField->runMBO();
        crossField->computeSingularities();

        // create tracer to get exact singularity locations
        SeparatrixTrace tracer(crossField, true);

        // collect theta and radius for up to 4 singularities
        double theta[4] = {0.0, 0.0, 0.0, 0.0};
        double radius[4] = {0.0, 0.0, 0.0, 0.0};
        int nSing = std::min(static_cast<int>(tracer.singularities.size()), 4);
        if (tracer.singularities.size() > 4) {
            std::cout << "ERROR: More than 4 singularities found, only recording the first 4\n";
            return 1;
        }
        for (int i = 0; i < nSing; ++i) {
            double x = tracer.singularities[i].coordinates[0];
            double y = tracer.singularities[i].coordinates[1];
            radius[i] = std::sqrt(x * x + y * y);
            theta[i] = std::atan2(y, x) * 180.0 / M_PI; // convert to degrees
            // convert theta to [0, 360)
            if (theta[i] < 0) {
                theta[i] += 360.0;
            }
        }

        // order theta and radius by theta for consistency
        for (int i = 0; i < nSing - 1; ++i) {
            for (int j = i + 1; j < nSing; ++j) {
                if (theta[j] < theta[i]) {
                    std::swap(theta[i], theta[j]);
                    std::swap(radius[i], radius[j]);
                }
            }
         }

        // write CSV row
        out << angle;
        for (int i = 0; i < 4; ++i) {
            out << "," << theta[i] << "," << radius[i];
        }
        out << "\n";

        // write VTK files for this frame
        std::ostringstream tag;
        tag << std::setw(4) << std::setfill('0') << frameIdx;
        std::string meshFile = "mesh_" + tag.str() + ".vtu";
        std::string singFile = "sing_" + tag.str() + ".vtp";
        std::string vtmFile  = "frame_" + tag.str() + ".vtm";

        std::cout << "  Writing VTK files: " << meshFile << ", " << singFile << ", " << vtmFile << "\n";
        writeMeshVTU(vtkDir + "/" + meshFile, *mesh);
        writeSingularitiesVTP(vtkDir + "/" + singFile, tracer.singularities);
        writeVTM(vtkDir + "/" + vtmFile, meshFile, singFile);

        pvdEntries.push_back({angle, vtkDir + "/" + vtmFile});

        // print to console
        std::cout << "  " << nSing << " singularities found\n";
    }

    out.close();

    // write PVD collection file (open in ParaView to play as movie)
    {
        std::ofstream pvd("circle_topologies.pvd");
        pvd << "<?xml version=\"1.0\"?>\n";
        pvd << "<VTKFile type=\"Collection\" version=\"1.0\">\n";
        pvd << "  <Collection>\n";
        for (const auto &entry : pvdEntries) {
            pvd << "    <DataSet timestep=\"" << entry.first
                << "\" file=\"" << entry.second << "\"/>\n";
        }
        pvd << "  </Collection>\n";
        pvd << "</VTKFile>\n";
        pvd.close();
    }

    std::cout << "CSV written to circle_topologies.csv\n";
    std::cout << "VTK files written to " << vtkDir << "/\n";
    std::cout << "Open circle_topologies.pvd in ParaView to view as movie\n";

    return 0;
}