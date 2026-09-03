// Utility to load a mesh, run the p=0 dual-mesh MBO cross-field solver,
// and print singularities and per-triangle field values.
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "dualMBO/DualMBO.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "Parameterization/CutMesh.hxx"
#include "Parameterization/UVGParam.hxx"

static const int MBO_MAX_STEPS = 500;

int main(int argc, char **argv) {
    // --no-pin turns off the disk-center Dirichlet pin, which is the state of
    // the solver before it existed: on a disk the four cones then land at an
    // arbitrary rotation instead of at 45 + k*90 degrees. It is here so the
    // two can be compared on the same mesh.
    bool pinDiskCenters = true;
    // On a multi-material mesh the interfaces are Dirichlet data for the field
    // in the way dS is, and a field that has not been told about them runs
    // straight through them and carries none of the cones the regions either
    // side of it need. MERIDIAN, TORSION and the viewer all set them; this
    // switch is how TestDualMBO gets at the same field, and it is what makes a
    // model whose disks are inclusions rather than the domain itself -- e.g.
    // multimat/bubbles -- show any cones inside those disks at all.
    bool alignInterfaces = false;
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--no-pin") pinDiskCenters = false;
        else if (a == "--align-interfaces") alignInterfaces = true;
        else pos.push_back(a);
    }

    if (pos.empty()) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [gamma] [max_steps] [--no-pin]\n";
        std::cerr << "  gamma     : edge penalty parameter (default 20.0)\n";
        std::cerr << "  max_steps : max MBO iterations   (default " << MBO_MAX_STEPS << ")\n";
        std::cerr << "  --no-pin  : do not pin disk centers (leaves the rotation free)\n";
        std::cerr << "  --align-interfaces : make material interfaces Dirichlet data too\n";
        return 1;
    }

    std::string path = pos[0];
    double gamma    = (pos.size() >= 2) ? std::stod(pos[1]) : 10.0;
    int    maxSteps = (pos.size() >= 3) ? std::stoi(pos[2]) : MBO_MAX_STEPS;

    // ------------------------------------------------------------------
    // Load mesh
    // ------------------------------------------------------------------
    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    // ------------------------------------------------------------------
    // Initialize DualMBO solver (assembles M, K, b and factorises A = M + tau*K)
    // ------------------------------------------------------------------
    DualMBO dualMBO(mesh, maxSteps, gamma);
    dualMBO.setPinDiskCenters(pinDiskCenters);

    std::unique_ptr<Interfaces> interfaces;
    if (alignInterfaces) {
        interfaces = std::make_unique<Interfaces>(mesh);
        if (interfaces->multiMaterial()) {
            dualMBO.setAlignedInteriorEdges(interfaces->interfaceEdges());
            std::cout << "Aligned to " << interfaces->interfaceEdges().size()
                      << " interface edge(s).\n";
        }
    }

    dualMBO.initialize();

    // ------------------------------------------------------------------
    // MBO iteration loop
    // ------------------------------------------------------------------
    int  stepCount = 0;
    bool converged = false;
    double ntris   = static_cast<double>(mesh->triangles.size());

    for (int i = 0; i < maxSteps; ++i) {
        dualMBO.step();
        ++stepCount;
        if (dualMBO.error < 2.0 * ntris * 1e-5) {
            converged = true;
            break;
        }
    }

    // ------------------------------------------------------------------
    // Detect singularities
    // ------------------------------------------------------------------
    dualMBO.computeSingularities();

    // ------------------------------------------------------------------
    // Disk report
    //
    // Every component whose boundary came out a circle, the triangle its
    // center pin landed on, and where the cones inside it ended up. The
    // interesting column is the last: the cone's polar angle about the disk
    // center reduced mod 90 degrees. A disk is rotationally symmetric, so
    // without the center pin that number is arbitrary and differs between
    // disks; with u pinned to 1 at the center the four +1/4 cones sit at
    // 45 + k*90 degrees, so every one of them should read close to 45.
    // ------------------------------------------------------------------
    {
        const auto &pinned = dualMBO.getDiskCenterTriangles();
        int nDisks = 0;
        for (const auto &c : mesh->materialComponents) if (c.circle.isCircle) ++nDisks;

        std::cout << "\nDisks: " << nDisks << " circular component(s), "
                  << pinned.size() << " center triangle(s) pinned to (1,0)\n";

        // vertex -> which component(s) it touches, for placing a singularity
        for (std::size_t ci = 0; ci < mesh->materialComponents.size(); ++ci) {
            const auto &comp = mesh->materialComponents[ci];
            if (!comp.circle.isCircle) continue;

            std::cout << "  disk (mat " << comp.matId << ") center ("
                      << std::fixed << std::setprecision(4)
                      << comp.circle.center[0] << ", " << comp.circle.center[1]
                      << ") r=" << comp.circle.radius
                      << "  centerTri=" << comp.centerTriangle << "\n";

            double sumIdx = 0.0;
            int nCones = 0;
            for (const auto &[v, idx] : dualMBO.singularVertices) {
                // a cone belongs to this disk if any incident triangle is in it
                const auto &vt = mesh->vertexTriangles;
                bool inside = false;
                for (int k = vt.rowPtr[v]; k < vt.rowPtr[v + 1]; ++k) {
                    if (mesh->triangleComponent[vt.colIdx[k]] == static_cast<int>(ci)) { inside = true; break; }
                }
                if (!inside) continue;

                Point d = mesh->vertices[v] - comp.circle.center;
                double ang = std::atan2(d[1], d[0]) * 180.0 / M_PI;
                if (ang < 0.0) ang += 360.0;
                double mod90 = std::fmod(ang, 90.0);
                double rfrac = normP(d) / comp.circle.radius;

                std::cout << "    cone v=" << v << " index=" << std::setprecision(3) << idx
                          << "  r/R=" << std::setprecision(3) << rfrac
                          << "  angle=" << std::setprecision(2) << ang
                          << " deg  (mod 90 = " << mod90 << ")\n";
                sumIdx += idx;
                ++nCones;
            }
            std::cout << "    -> " << nCones << " cone(s), index sum " << std::setprecision(3) << sumIdx << "\n";
        }
        std::cout << std::defaultfloat << "\n";
    }

    // ------------------------------------------------------------------
    // Cut mesh from DualMBO field
    // ------------------------------------------------------------------
    CutMesh cm(dualMBO);
    const auto &rep = cm.sanityCheck();

    if (rep.looksLikeDisk && rep.allSingularitiesOnBoundary)
        std::cout << "\033[32m[PASS]\033[0m Cut mesh is a disk with singularities on boundary.\n";
    else
        std::cout << "\033[31m[FAIL]\033[0m Cut mesh sanity check failed.\n";

    // ------------------------------------------------------------------
    // Global UV parameterization
    // ------------------------------------------------------------------
    try {
        UVGParam uvp(cm);

        int flips = uvp.numFlippedTriangles();
        if (flips == 0)
            std::cout << "\033[32m[PASS]\033[0m No flipped triangles in UV space.\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m " << flips << " flipped triangle(s) in UV space.\n";
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m UVGParam failed: " << e.what() << "\n";
    }

    return 0;
}