cmake_minimum_required(VERSION 3.15)
project(CrossGen LANGUAGES CXX)

# Build options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Eigen (header-only) - use local submodule
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen")
add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}"
)
message(STATUS "Using local Eigen3 from: ${EIGEN3_INCLUDE_DIR}")

# OpenMP
#
# Used by the TMOP smoother (src/mesh/TMOP.cxx) and by the Shape-DNA module
# (src/ShapeDNA: assembly, the nested-dissection Cholesky, the Lanczos
# eigensolver). The smoother sweeps in colour order and every ShapeDNA loop
# writes disjoint outputs or reduces in a fixed order, so neither answer depends
# on the thread count, and a build without OpenMP is a correct build, only a
# slower one -- every pragma is compiled out through _OPENMP.
#
# AppleClang is the awkward case: it understands the pragmas but ships no
# runtime, so plain find_package(OpenMP) fails on a Mac that has Homebrew's
# libomp installed and working. Fall back to wiring libomp up by hand there.
option(CROSSGEN_ENABLE_OPENMP "Use OpenMP to parallelise TMOP and the Shape-DNA solver" ON)

set(CROSSGEN_OPENMP_TARGET "")
if (CROSSGEN_ENABLE_OPENMP)
    find_package(OpenMP QUIET COMPONENTS CXX)
    if (OpenMP_CXX_FOUND)
        set(CROSSGEN_OPENMP_TARGET OpenMP::OpenMP_CXX)
        message(STATUS "OpenMP found: ${OpenMP_CXX_FLAGS}")
    elseif (APPLE)
        execute_process(COMMAND brew --prefix libomp
                        OUTPUT_VARIABLE LIBOMP_PREFIX
                        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        if (NOT LIBOMP_PREFIX)
            set(LIBOMP_PREFIX "/opt/homebrew/opt/libomp")
        endif()
        if (EXISTS "${LIBOMP_PREFIX}/include/omp.h")
            add_library(CrossGenOpenMP INTERFACE)
            target_compile_options(CrossGenOpenMP INTERFACE -Xpreprocessor -fopenmp)
            target_include_directories(CrossGenOpenMP INTERFACE "${LIBOMP_PREFIX}/include")
            if (EXISTS "${LIBOMP_PREFIX}/lib/libomp.dylib")
                target_link_libraries(CrossGenOpenMP INTERFACE "${LIBOMP_PREFIX}/lib/libomp.dylib")
            else()
                target_link_libraries(CrossGenOpenMP INTERFACE "${LIBOMP_PREFIX}/lib/libomp.a")
            endif()
            set(CROSSGEN_OPENMP_TARGET CrossGenOpenMP)
            message(STATUS "OpenMP via Homebrew libomp at ${LIBOMP_PREFIX}")
        endif()
    endif()
endif()
if (NOT CROSSGEN_OPENMP_TARGET)
    message(STATUS "No OpenMP: TMOP and ShapeDNA will run single-threaded")
endif()

# OpenCASCADE: the geometry kernel
#
# Every curve, surface and boundary representation in the repository is an
# OpenCASCADE object, and src/geom is the only code allowed to know that. It is
# built as its own library with OpenCASCADE linked PRIVATE, so the kernel's
# include directories reach src/geom's sources and nothing else: a file
# anywhere else that includes an OpenCASCADE header does not compile. ctest's
# Geom_OpenCascadeConfinedToGeom checks the same rule from the source side.
# (Homebrew: `brew install opencascade`; elsewhere, point OpenCASCADE_DIR at
# the directory holding OpenCASCADEConfig.cmake.)
find_package(OpenCASCADE REQUIRED)
message(STATUS "OpenCASCADE ${OpenCASCADE_VERSION} from ${OpenCASCADE_INSTALL_PREFIX}")

add_library(CrossGenGeom SHARED
    src/geom/Vec.hxx
    src/geom/BSpline.hxx
    src/geom/BSpline.cxx
    src/geom/Polyline.hxx
    src/geom/Polyline.cxx
    src/geom/ArcLength.hxx
    src/geom/ArcLength.cxx
    src/geom/Fitting.hxx
    src/geom/Fitting.cxx
    src/geom/Coons.hxx
    src/geom/Coons.cxx
    src/geom/Topology.hxx
    src/geom/Topology.cxx
    src/geom/detail/Occ.hxx
)
target_include_directories(CrossGenGeom PUBLIC src)
target_link_libraries(CrossGenGeom PRIVATE
    TKernel TKMath TKG2d TKG3d TKGeomBase TKGeomAlgo TKBRep TKTopAlgo
    TKXSBase TKDE TKDESTEP)

option(BUILD_MESH2GMSH "Build the Mesh2Dgmsh utility" ON)

if (BUILD_MESH2GMSH)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
    find_package(Gmsh REQUIRED)
    message(STATUS "Gmsh found: ${GMSH_LIBRARY}")
    add_executable(Mesh2Dgmsh src/utils/Mesh2Dgmsh.cxx)
    target_link_libraries(Mesh2Dgmsh PRIVATE Gmsh::gmsh)

    add_executable(Mesh3Dto2Dfacesgmsh src/utils/Mesh3Dto2Dfacesgmsh.cxx)
    target_link_libraries(Mesh3Dto2Dfacesgmsh PRIVATE Gmsh::gmsh)
else()
    message(STATUS "Skipping Mesh2Dgmsh utility")
endif()

add_library(PolyVector SHARED
    src/mesh/VertexTriangleCSR.hxx
    src/mesh/VertexTriangleCSR.cxx
    src/mesh/Mesh.cxx
    src/mesh/Mesh.hxx
    src/mesh/QuadMesh.cxx
    src/mesh/QuadMesh.hxx
    src/mesh/TMOP.cxx
    src/mesh/TMOP.hxx
    src/mesh/BlockDecomposition.hxx
    src/mesh/BlockDecomposition.cxx
    src/mesh/BlockQuadMesh.hxx
    src/mesh/BlockQuadMesh.cxx
    src/ShapeDNA/Parallel.hxx
    src/ShapeDNA/SparseMatrix.hxx
    src/ShapeDNA/SparseMatrix.cxx
    src/ShapeDNA/DenseEigen.hxx
    src/ShapeDNA/DenseEigen.cxx
    src/ShapeDNA/SparseCholesky.hxx
    src/ShapeDNA/SparseCholesky.cxx
    src/ShapeDNA/Lanczos.hxx
    src/ShapeDNA/Lanczos.cxx
    src/ShapeDNA/LagrangeFEM.hxx
    src/ShapeDNA/LagrangeFEM.cxx
    src/ShapeDNA/ShapeDNA.hxx
    src/ShapeDNA/ShapeDNA.cxx
    src/polyvector/PolyVectors.hxx
    src/polyvector/PolyVectors.cxx
    src/dualmbo/DualMBO.hxx
    src/dualmbo/DualMBO.cxx
    src/Parameterization/CutMesh.hxx
    src/Parameterization/CutMesh.cxx
    src/Parameterization/HarmonicCut.hxx
    src/Parameterization/HarmonicCut.cxx
    src/Parameterization/UVGParam.hxx
    src/Parameterization/UVGParam.cxx
    src/crossfield/CrossField.hxx
    src/crossfield/CrossField.cxx
    src/OASIS/OASIS.hxx
    src/OASIS/OASIS.cxx
    src/UMBER/UMBER.hxx
    src/UMBER/UMBER.cxx
    src/MERIDIAN/MERIDIAN.hxx
    src/MERIDIAN/MERIDIAN.cxx
    src/MERIDIAN/Interfaces.hxx
    src/MERIDIAN/Interfaces.cxx
    src/MERIDIAN/ConeSingularities.hxx
    src/MERIDIAN/ConeSingularities.cxx
    src/MERIDIAN/ConeCut.hxx
    src/MERIDIAN/ConeCut.cxx
    src/MERIDIAN/RicciFlow.hxx
    src/MERIDIAN/RicciFlow.cxx
    src/MERIDIAN/Immersion.hxx
    src/MERIDIAN/Immersion.cxx
    src/MERIDIAN/SubdomainLabels.hxx
    src/MERIDIAN/SubdomainLabels.cxx
    src/MERIDIAN/LayoutEnergy.hxx
    src/MERIDIAN/LayoutEnergy.cxx
    src/MERIDIAN/Separatrices.hxx
    src/MERIDIAN/Separatrices.cxx
    src/MERIDIAN/Arrangement.hxx
    src/MERIDIAN/Arrangement.cxx
    src/MERIDIAN/DiskTemplate.hxx
    src/MERIDIAN/DiskTemplate.cxx
    src/MERIDIAN/SplineFit.hxx
    src/MERIDIAN/SplineFit.cxx
    src/MERIDIAN/QuadMesh.hxx
    src/MERIDIAN/QuadMesh.cxx
    src/TORSION/ConeMetric.hxx
    src/TORSION/ConeMetric.cxx
    src/TORSION/FieldFrames.hxx
    src/TORSION/FieldFrames.cxx
    src/TORSION/TutteEmbedding.hxx
    src/TORSION/TutteEmbedding.cxx
    src/TORSION/FieldIntegration.hxx
    src/TORSION/FieldIntegration.cxx
    src/TORSION/TORSION.hxx
    src/TORSION/TORSION.cxx
    src/ATLAS/SquareTransport.hxx
    src/ATLAS/PlanarDomain.hxx
    src/ATLAS/PlanarDomain.cxx
    src/ATLAS/SquareCarrier.hxx
    src/ATLAS/SquareCarrier.cxx
    src/ATLAS/CavityFill.hxx
    src/ATLAS/CavityFill.cxx
    src/ATLAS/ExplicitTemplates.hxx
    src/ATLAS/ExplicitTemplates.cxx
    src/ATLAS/RectangleCertifier.hxx
    src/ATLAS/RectangleCertifier.cxx
    src/ATLAS/BlockCover.hxx
    src/ATLAS/BlockCover.cxx
    src/ATLAS/CavityRewrite.hxx
    src/ATLAS/CavityRewrite.cxx
    src/ATLAS/CoarseDomain.hxx
    src/ATLAS/CoarseDomain.cxx
    src/ATLAS/ReferenceField.hxx
    src/ATLAS/ReferenceField.cxx
    src/ATLAS/Realisation.hxx
    src/ATLAS/Realisation.cxx
    src/ATLAS/BlockMesh.hxx
    src/ATLAS/BlockMesh.cxx
    src/ATLAS/ATLAS.hxx
    src/ATLAS/ATLAS.cxx
    src/UMBER/Polysquare.hxx
    src/UMBER/Polysquare.cxx
    src/UMBER/MotorcycleGraph.hxx
    src/UMBER/MotorcycleGraph.cxx
    src/UMBER/BlockLayout.hxx
    src/UMBER/BlockLayout.cxx
    src/ZIPLINE/FieldTracer.hxx
    src/ZIPLINE/FieldTracer.cxx
    src/ZIPLINE/SeparatrixTrace.hxx
    src/ZIPLINE/SeparatrixTrace.cxx
    src/ZIPLINE/QuadLayout.hxx
    src/ZIPLINE/QuadLayout.cxx
    src/ZIPLINE/PartitionSimplify.hxx
    src/ZIPLINE/PartitionSimplify.cxx
    src/ZIPLINE/LayoutBlocks.hxx
    src/ZIPLINE/LayoutBlocks.cxx
    src/ZIPLINE/StemExtension.hxx
    src/ZIPLINE/StemExtension.cxx
    src/ZIPLINE/TriangleLocator.hxx
    src/ZIPLINE/ZIPLINE.hxx
    src/ZIPLINE/ZIPLINE.cxx
    src/triangle/triangle.cpp
    src/triangle/triangle.h
    src/triangle/TriangleMesher.hpp
    src/medialaxis/MedialAxis.hxx
    src/medialaxis/MedialAxis.cxx
    src/medialaxis/MedialAxisMap.hxx
    src/medialaxis/MedialAxisMap.cxx
    src/medialaxis/MedialAxisTMesh.hxx
    src/medialaxis/MedialAxisTMesh.cxx
    src/quantization/QuantTMesh.hxx
    src/quantization/QuantTMesh.cxx
    src/quantization/TMeshQuantizer.hxx
    src/quantization/TMeshQuantizer.cxx
    src/quantization/QuantTMeshConvert.hxx
    src/quantization/QuantTMeshConvert.cxx
    src/quantization/TMeshContract.hxx
    src/quantization/TMeshContract.cxx
)
target_include_directories(PolyVector PUBLIC src)
target_include_directories(PolyVector PUBLIC src/triangle)
target_compile_definitions(PolyVector PRIVATE TRILIBRARY ANSI_DECLARATORS)
target_link_libraries(PolyVector PUBLIC Eigen3::Eigen CrossGenGeom)
if (CROSSGEN_OPENMP_TARGET)
    target_link_libraries(PolyVector PUBLIC ${CROSSGEN_OPENMP_TARGET})
endif()

# Small utility to load and print mesh data
add_executable(LoadMesh src/utils/LoadMesh.cxx)
target_link_libraries(LoadMesh PRIVATE PolyVector)

# ZIPLINE, the Viertel et al. 2019 tracing method (src/ZIPLINE): separatrices,
# the T-layout they cut out, its Sec. 4 simplification, blocks and a mesh
add_executable(TestZIPLINE src/utils/TestZIPLINE.cxx)
target_link_libraries(TestZIPLINE PRIVATE PolyVector)

# Utility to cut a mesh into a disk and run sanity checks
add_executable(TestCutMesh src/utils/TestCutMesh.cxx)
target_link_libraries(TestCutMesh PRIVATE PolyVector)

# Test medial axis construction and phi-inverse mapping
add_executable(TestMedialAxis src/utils/TestMedialAxis.cxx)
target_link_libraries(TestMedialAxis PRIVATE PolyVector)

# Test the QGP T-mesh quantizer and its converters
add_executable(TestQuantization src/utils/TestQuantization.cxx)
target_link_libraries(TestQuantization PRIVATE PolyVector)

# Generate circle/ellipse topologies
add_executable(GenerateCircleTopologies src/utils/GenerateCircleTopologies.cxx)
target_link_libraries(GenerateCircleTopologies PRIVATE PolyVector)

# Test OASIS quasi-eigenfunction on a meshed rectangle
add_executable(TestOasis src/utils/TestOasis.cxx)
target_link_libraries(TestOasis PRIVATE PolyVector)

# Test p=0 dual-mesh MBO cross-field solver
add_executable(TestDualMBO src/utils/TestDualMBO.cxx)
target_link_libraries(TestDualMBO PRIVATE PolyVector)

# Smooth a quadrilateral mesh with the in-house TMOP optimizer of
# src/mesh/TMOP.hxx, and self-test it together with the node-mobility layer of
# mesh::QuadMesh it stands on. Takes a *quad* .obj -- what TestMERIDIAN --mesh
# writes; --selftest builds its own meshes and needs no files.
add_executable(TestTMOP src/utils/TestTMOP.cxx)
target_link_libraries(TestTMOP PRIVATE PolyVector)

# Self-test the geometry layer of src/geom -- B-spline curves and surfaces,
# fitting, interpolation, knot insertion, arc length, projection, Coons blends
# and the B-rep of vertices, edges and faces -- on shapes whose answers are
# known in closed form. Needs no files, and links src/geom alone.
add_executable(TestGeom src/utils/TestGeom.cxx)
target_link_libraries(TestGeom PRIVATE CrossGenGeom)

# Test the UMBER frame field optimization (Wang et al. 2022, Secs. 4.1-4.2)
add_executable(TestUMBER src/utils/TestUMBER.cxx)
target_link_libraries(TestUMBER PRIVATE PolyVector)
# The driver writes its VTK output into the working directory, so it runs in
# the build tree rather than beside the meshes.
set_target_properties(TestUMBER PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Test MERIDIAN Stages 1-9 (Shepherd, Gu and Hughes 2022, Secs. 3.1-3.3, 4 and
# 5): cone singularities, cutting graph, Gauss-Bonnet, discrete Ricci flow, the
# metric immersion, the layout energies, the separatrices of the result, the
# arrangement they cut S into, and the bicubic patches fitted to it
add_executable(TestMERIDIAN src/utils/TestMERIDIAN.cxx)
target_link_libraries(TestMERIDIAN PRIVATE PolyVector)

# Test TORSION, Pipeline B of docs/cf_flow_pipeline.md: the same Stages 1-10 as
# MERIDIAN with psi_0 built by integrating the DualMBO cross field instead of by
# unfolding a Ricci metric. --selftest exercises the four new pieces -- the
# combing and its matchings, the Tutte embedding, the constrained integration
# and Reference::Field -- on meshes whose answers are known in closed form.
add_executable(TestTORSION src/utils/TestTORSION.cxx)
target_link_libraries(TestTORSION PRIVATE PolyVector)

# The Shape-DNA of a planar mesh (src/ShapeDNA, docs/shape_dna.md): Lagrange
# P1-P3 elements, a nested-dissection Cholesky and a shift-invert block Lanczos
# eigensolver, all in-house and OpenMP-parallel. Given meshes, prints each one's
# DNA and their pairwise distances; --selftest checks the module against the
# rectangle and disk spectra and needs no files.
add_executable(TestShapeDNA src/utils/TestShapeDNA.cxx)
target_link_libraries(TestShapeDNA PRIVATE PolyVector)

# Test ATLAS, the square-transport blocking of
# docs/square_transport_2d_theory_and_implementation.md, Stages 1-6: the
# three-quad carrier and its integer transports, whole-region templates,
# certified rectangles, a conforming cover, and cavity rewrites. --selftest
# runs the spec's regression geometries whose answers are known in closed form.
add_executable(TestATLAS src/utils/TestATLAS.cxx)
target_link_libraries(TestATLAS PRIVATE PolyVector)

# copy all meshes to build directory for easy access, keeping the material
# subdirectory layout (singlemat/, multimat/, ...) intact
file(GLOB_RECURSE MESH_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/data/meshes"
     "${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/*")
foreach(MESH_REL ${MESH_FILES})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/${MESH_REL}
                   ${CMAKE_CURRENT_BINARY_DIR}/data/meshes/${MESH_REL} COPYONLY)
endforeach()

# Custom target to run TestCutMesh on all meshes
file(GLOB_RECURSE MESH_OBJ_FILES "${CMAKE_CURRENT_BINARY_DIR}/data/meshes/*.obj")
set(CUT_COMMANDS "")
foreach(MESH_FILE ${MESH_OBJ_FILES})
    get_filename_component(MESH_NAME ${MESH_FILE} NAME)
    list(APPEND CUT_COMMANDS 
        COMMAND ${CMAKE_COMMAND} -E echo "Processing ${MESH_NAME}..."
        COMMAND $<TARGET_FILE:TestCutMesh> ${MESH_FILE}
    )
endforeach()

add_custom_target(cut
    ${CUT_COMMANDS}
    DEPENDS TestCutMesh
    COMMENT "Running TestCutMesh on all meshes in data/meshes"
    VERBATIM
)

# Custom target to run TestDualMBO on all meshes
set(DualMBO_COMMANDS "")
foreach(MESH_FILE ${MESH_OBJ_FILES})
    get_filename_component(MESH_NAME ${MESH_FILE} NAME)
    list(APPEND DualMBO_COMMANDS 
        COMMAND ${CMAKE_COMMAND} -E echo "Processing ${MESH_NAME}..."
        COMMAND $<TARGET_FILE:TestDualMBO> ${MESH_FILE} || ${CMAKE_COMMAND} -E true
    )
endforeach()

add_custom_target(dualmbo
    ${DualMBO_COMMANDS}
    DEPENDS TestDualMBO
    COMMENT "Running TestDualMBO on all meshes in data/meshes"
    VERBATIM
)

# Custom target to run TestMERIDIAN on all meshes
set(MERIDIAN_COMMANDS "")
foreach(MESH_FILE ${MESH_OBJ_FILES})
    get_filename_component(MESH_NAME ${MESH_FILE} NAME)
    list(APPEND MERIDIAN_COMMANDS
        COMMAND ${CMAKE_COMMAND} -E echo "Processing ${MESH_NAME}..."
        COMMAND $<TARGET_FILE:TestMERIDIAN> ${MESH_FILE} || ${CMAKE_COMMAND} -E true
    )
endforeach()

add_custom_target(meridian
    ${MERIDIAN_COMMANDS}
    DEPENDS TestMERIDIAN
    COMMENT "Running TestMERIDIAN on all meshes in data/meshes"
    VERBATIM
)

# Custom target to run TestTORSION on all meshes
set(TORSION_COMMANDS "")
foreach(MESH_FILE ${MESH_OBJ_FILES})
    get_filename_component(MESH_NAME ${MESH_FILE} NAME)
    list(APPEND TORSION_COMMANDS
        COMMAND ${CMAKE_COMMAND} -E echo "Processing ${MESH_NAME}..."
        COMMAND $<TARGET_FILE:TestTORSION> ${MESH_FILE} || ${CMAKE_COMMAND} -E true
    )
endforeach()

add_custom_target(torsion
    ${TORSION_COMMANDS}
    DEPENDS TestTORSION
    COMMENT "Running TestTORSION on all meshes in data/meshes"
    VERBATIM
)

# Enable testing
enable_testing()

# Test executable for TriangleMesher
add_executable(TestTriangleMesher src/triangle/tests/TestTriangleMesher.cxx)
target_link_libraries(TestTriangleMesher PRIVATE PolyVector)

# Add individual CTest tests
add_test(NAME TriangleMesher_Rectangle COMMAND TestTriangleMesher 1)
add_test(NAME TriangleMesher_Circle COMMAND TestTriangleMesher 2)
add_test(NAME TriangleMesher_CircleWithEllipseHole COMMAND TestTriangleMesher 3)
add_test(NAME TriangleMesher_CircleWithEllipseRegion COMMAND TestTriangleMesher 4)
add_test(NAME TriangleMesher_CircleWithTwoEllipseHoles COMMAND TestTriangleMesher 5)
add_test(NAME TriangleMesher_ConcentricShells COMMAND TestTriangleMesher 6)


# End-to-end: trace every model in data/meshes and check each layout is a plane
# graph whose faces tile the model exactly with nothing left dangling, that the
# blocks read off it (ZIPLINE/LayoutBlocks) are a sound decomposition -- every
# side running between its block's corners, every Coons patch on its corners --
# and that the mesh on them is conforming. How many of the components come out
# four-sided, and how much of the model the blocks cover, are printed but not
# asserted: that is what the method achieves on a model, not a defect.
add_test(NAME ZIPLINE_LayoutsOnAllMeshes
         COMMAND TestZIPLINE --no-vtu --mesh 0.05 ${MESH_OBJ_FILES})

# docs/viertel_2019.md Sec. 11 on closed-form fields: the index, the ports and
# the hyperbola of a singularity, straight streamlines in a constant field, and
# the square, the L-shape and a T-junction taken through to blocks and a mesh.
add_test(NAME ZIPLINE_SelfTest COMMAND TestZIPLINE --selftest)

# Discrete Ricci flow (Shepherd et al. Sec. 3.2.1) against prescribed cone sets
# with far more curvature to move than anything in data/meshes asks for, which
# is what puts the weighted-Delaunay flipping under load. Checks the cone angle
# sums directly and the Hessian against a finite difference.
add_test(NAME MERIDIAN_RicciFlowSelfTest COMMAND TestMERIDIAN --selftest)

# The geometry layer under both pipelines' Stage 9 (src/geom), against closed
# forms: partition of unity, knot insertion and reversal leaving the curve
# unchanged, fits reproducing what they can represent, a Coons net on four
# straight sides reproducing the bilinear map, and two faces sharing an edge
# making a valid shape with that edge counted once.
add_test(NAME Geom_SplineSelfTest COMMAND TestGeom)

# The rule the CrossGenGeom target enforces at compile time, checked on the
# sources: no file outside src/geom includes an OpenCASCADE header.
add_test(NAME Geom_OpenCascadeConfinedToGeom
         COMMAND ${CMAKE_COMMAND} -DROOT=${CMAKE_CURRENT_SOURCE_DIR}
                 -DOCC_INCLUDE=${OpenCASCADE_INCLUDE_DIR}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CheckOpenCascadeConfinement.cmake)

# End-to-end: Stages 1 to 7 on a handful of models that between them exercise
# every part of the later stages -- a plain disk with one interior cone
# (geom001), a model whose cutting graph carries several different holonomies
# (geom003), one with voids and a dozen connectivity constraints (geom012), and
# a two-material model, where the interface between the physical surfaces is a
# feature chain and E3 has something to hold (multimat geom001). Each asserts
# Q1 to Q5 on the map that comes out and traces its separatrices; the whole
# corpus is `make meridian`.
#
# The exit code follows the layout, not the tracing. A separatrix that runs on
# past the step cap is Q5 failing on a curve E5 was never given a constraint
# for -- Remark 3.1's case, and one whose remedy is a Gamma_topo constraint
# rather than a fix to Stage 7 -- so it is reported as a warning. geom012 has
# four of them and passes.
foreach(MERIDIAN_CASE singlemat/geom001 singlemat/geom003 singlemat/geom012 multimat/geom001)
    string(REPLACE "/" "_" MERIDIAN_CASE_NAME ${MERIDIAN_CASE})
    add_test(NAME MERIDIAN_Layout_${MERIDIAN_CASE_NAME}
             COMMAND TestMERIDIAN ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/${MERIDIAN_CASE}.obj)
endforeach()

# Stages 0c and 11, the circular-inclusion templates, on the model they exist
# for: ten disks in a matrix block. Without them the layout has eighty cones it
# cannot place and leaves 57% of the model unmeshed; with them the excised
# problem reaches Definition 2.1 and every inclusion comes back as an O-grid.
# See DiskTemplate.
#
# --near-miss is not decoration. At the default 0.15 Gamma_topo over-seeds the
# matrix and four of the ten rims come back without a closed loop for the
# template to attach to; 0.04 is the tolerance measured in
# SubdomainLabels::Options::nearMissTolerance and it is what makes this model
# converge. The exit code is the whole pipeline's, templates included.
add_test(NAME MERIDIAN_Templates_multimat_bubbles
         COMMAND TestMERIDIAN ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/multimat/bubbles.obj
                 --disk-templates --near-miss 0.04)
set_tests_properties(MERIDIAN_Templates_multimat_bubbles PROPERTIES TIMEOUT 2400)

# The four pieces Pipeline B adds -- the combing and its matchings, the
# constrained integration, the Tutte embedding, and the field reference -- on
# inputs whose answers are known in closed form before the code runs. A square
# whose boundary is axis aligned has a *constant* cross field: every matching is
# zero, the integration's answer is exactly x -> x/h, and Reference::Field must
# report det J = 1 on every triangle. A circle has none of that and is there for
# the branch selection and the seam transitions to have something to do.
# Case 3 adds Sec. 4's cone metric on both: the square's four corners already
# carry the pi/2 their index asks for, so u is constant there and any movement
# is movement for no reason; the disk has 2 pi of boundary turning to move and
# what is known of it is the target, K = Kbar to rounding on a metric that is
# still realisable.
# UMBER end to end: the frame field of Sec. 4.2, the polysquare of Sec. 4.3,
# the motorcycle graph of Sec. 5, the chord collapse, and then the two stages
# every method now shares -- the block decomposition (mesh/BlockDecomposition)
# and the quadrilateral mesh on it (mesh/BlockQuadMesh).
#
# What the exit code asserts is the structural part and only that: every side
# of every block resolves to a macro edge, every interface macro edge separates
# two materials, no block sits in two materials, and the mesh is conforming.
# The flips, the internal singularities and the boundary turn count are printed
# and not asserted, because whether a given shape has a polysquare at all is a
# fact about the shape.
#
# Two models: geom012, the heaviest single-material layout in the corpus (ten
# holes, 107 blocks), and multimat/geom001, the quarter disk in a square, which
# is the smallest domain on which the interfaces have to be features of the
# field, the deformation and the block structure alike.
add_test(NAME UMBER_Polysquare_singlemat_geom012
         COMMAND TestUMBER ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/singlemat/geom012.obj)
add_test(NAME UMBER_Polysquare_multimat_geom001
         COMMAND TestUMBER ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/multimat/geom001.obj)
set_tests_properties(UMBER_Polysquare_singlemat_geom012
                     UMBER_Polysquare_multimat_geom001 PROPERTIES TIMEOUT 2400)

add_test(NAME TORSION_SelfTest COMMAND TestTORSION --selftest)

# End-to-end on the same four models MERIDIAN_Layout uses, with psi_0 integrated
# from the cross field instead of unfolded from a Ricci metric. Asserts Q1-Q5 on
# the map that comes out, exactly as the MERIDIAN cases do; what is *not*
# asserted is the arrangement, because the two routes differ there for a reason
# that is structural rather than numerical -- see FieldFrames::Report::
# boundaryTurnMismatches.
foreach(TORSION_CASE singlemat/geom001 singlemat/geom003 singlemat/geom012 multimat/geom001)
    string(REPLACE "/" "_" TORSION_CASE_NAME ${TORSION_CASE})
    add_test(NAME TORSION_Layout_${TORSION_CASE_NAME}
             COMMAND TestTORSION ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/${TORSION_CASE}.obj)
    # Longer than the MERIDIAN cases need. The continuation of Stage 6 stops
    # early when every property of Definition 2.1 holds, and on a
    # field-integrated map it often does not get to: the boundary corners
    # FieldFrames predicts sit at vertices Stage 1 called regular, so Q2 has a
    # residual the penalties cannot spend away and the full schedule runs.
    set_tests_properties(TORSION_Layout_${TORSION_CASE_NAME} PROPERTIES TIMEOUT 2400)
endforeach()

# The same templates through Pipeline B. Stages 0c and 11 are shared code and
# the point of running them both ways is that the excised problem is the same
# problem: both routes reach a layout on it and both fill all ten inclusions.
add_test(NAME TORSION_Templates_multimat_bubbles
         COMMAND TestTORSION ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/multimat/bubbles.obj
                 --disk-templates --near-miss 0.04)
set_tests_properties(TORSION_Templates_multimat_bubbles PROPERTIES TIMEOUT 2400)

# Sec. 4 of docs/cf_flow_pipeline.md asks for the field reference to be tested
# standalone, from psi_R, before any integration exists to be blamed for its
# failures. This is that test, over all four references at once, and its
# verdicts are the finding it produced: Reference::Field is Reference::Euclidean
# to machine precision, because the frame is a scaled rotation and the symmetric
# Dirichlet energy cannot see a rotation on the right of J. geom003 is the model
# Sec. 10 names for it -- the disk whose cone fan the Euclidean reference
# shreds.
add_test(NAME TORSION_ReferenceMetrics
         COMMAND TestTORSION --ref-test
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/singlemat/geom003.obj)
set_tests_properties(TORSION_ReferenceMetrics PROPERTIES TIMEOUT 2400)

# ATLAS on inputs whose blocking is known before it runs (Sec. 13.3 of the
# spec): one triangle (Sec. 3.1's corner determinants exactly, three blocks), a
# rectangle (one block), a disk (a five-block O-grid), a half disk (a four-block
# half O-grid, one block at each corner), an annulus (four sectors), an L
# (three blocks), two materials (two blocks sharing the
# interface), the certificate refusing a periodic band and an L-tromino,
# Stage 6 alone lowering the valence defect of a raw carrier, a plate with two
# holes blocked by a coarse search and realised on the input, two meshes of
# one boundary giving blockings of about the same size, and the TFI mesh on the
# blocks (BlockMesh) where its counts and boundary are known.
add_test(NAME ATLAS_SelfTest COMMAND TestATLAS --selftest)

# The TFI mesh on the blocks and TMOP after it, as the viewer's ATLAS mode runs
# them; exit 9 if the smoothed mesh folds. geom014 is where TMOP sampled at the
# 2x2 Gauss points turned four corners of a valid TFI mesh over (worst scaled
# Jacobian 0.19 -> -0.81); sampled at the corners it reaches 0.46 with none.
add_test(NAME ATLAS_Mesh_singlemat_geom014
         COMMAND TestATLAS ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/singlemat/geom014.obj
                 --mesh 0.05 --mesh-tmop 1000)
set_tests_properties(ATLAS_Mesh_singlemat_geom014 PROPERTIES TIMEOUT 600)

# The corpus's half disk. Before the half O-grid it was a five-block O-grid with
# each 90-degree corner cut into two cells and a worst TFI element of 0.054; now
# four blocks and 0.70. Exit 8 above four blocks, 9 if the mesh does not validate.
add_test(NAME ATLAS_HalfOGrid_singlemat_geom002
         COMMAND TestATLAS ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/singlemat/geom002.obj
                 --max-blocks 4 --mesh 0.05)
set_tests_properties(ATLAS_HalfOGrid_singlemat_geom002 PROPERTIES TIMEOUT 600)

# The coarse searches on models no template covers whole -- holes, a slot, a
# gear -- with a ceiling on the block count well above what they reach now
# (13-20, 8-12, 55-70 on 2026-09-18) and far below the thousands the fine
# carrier alone gives. A regression in the search shows up as a failure here
# rather than as a still-valid fine fallback.
foreach(ATLAS_CASE "singlemat/geom007 30" "singlemat/geom021 20" "singlemat/geom009 100")
    string(REPLACE " " ";" ATLAS_ARGS ${ATLAS_CASE})
    list(GET ATLAS_ARGS 0 ATLAS_MODEL)
    list(GET ATLAS_ARGS 1 ATLAS_MAX)
    string(REPLACE "/" "_" ATLAS_CASE_NAME ${ATLAS_MODEL})
    add_test(NAME ATLAS_Coarse_${ATLAS_CASE_NAME}
             COMMAND TestATLAS ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/${ATLAS_MODEL}.obj --max-blocks ${ATLAS_MAX})
    set_tests_properties(ATLAS_Coarse_${ATLAS_CASE_NAME} PROPERTIES TIMEOUT 600)
endforeach()

# The same, on the multi-material models: until CoarseDomain carried the
# interface network these had no coarse search at all and fell back on the
# fine carrier, where Stage 6 stalls (geom001 29054 blocks against the 8 its
# layout wants, geom012 43406, the corpus 225617). With it, each is a handful
# to a few dozen. The ceilings are about three times what they reach now
# (12, 12, 9, 30 on 2026-09-19), so annealer seed noise cannot trip them but a
# lost coarse search will.
foreach(ATLAS_CASE "multimat/geom001 40" "multimat/geom003 40" "multimat/geom012 40" "multimat/geom013 90")
    string(REPLACE " " ";" ATLAS_ARGS ${ATLAS_CASE})
    list(GET ATLAS_ARGS 0 ATLAS_MODEL)
    list(GET ATLAS_ARGS 1 ATLAS_MAX)
    string(REPLACE "/" "_" ATLAS_CASE_NAME ${ATLAS_MODEL})
    add_test(NAME ATLAS_Coarse_${ATLAS_CASE_NAME}
             COMMAND TestATLAS ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/${ATLAS_MODEL}.obj --max-blocks ${ATLAS_MAX})
    set_tests_properties(ATLAS_Coarse_${ATLAS_CASE_NAME} PROPERTIES TIMEOUT 600)
endforeach()

# The reference cross field (docs/atlas_crossfield_guidance.md). geom025, a
# rounded plate with a V notch, is where ATLAS put three of the field's four
# interior +1/4 cones on the arcs as block corners: 4 blocks, each corner one
# cell spanning pi, TMOP worst 0.25-0.36 over three seeds. With the field and
# the sign-aware boundary defects the cones stay inside: 12 blocks, TMOP worst
# 0.82-0.85. Exit 8 above 16 blocks, 10 below 0.6.
add_test(NAME ATLAS_Field_singlemat_geom025
         COMMAND TestATLAS ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/singlemat/geom025.obj
                 --max-blocks 16 --mesh-frac 0.025 --mesh-tmop 100 --min-sj 0.6)
set_tests_properties(ATLAS_Field_singlemat_geom025 PROPERTIES TIMEOUT 600)

# End-to-end on a single- and a multi-material model. The exit code is the
# validity of the blocking, which the fine fallback also satisfies; how coarse
# it came out is printed, not asserted.
foreach(ATLAS_CASE singlemat/geom001 multimat/geom001)
    string(REPLACE "/" "_" ATLAS_CASE_NAME ${ATLAS_CASE})
    add_test(NAME ATLAS_Blocking_${ATLAS_CASE_NAME}
             COMMAND TestATLAS ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/${ATLAS_CASE}.obj)
    set_tests_properties(ATLAS_Blocking_${ATLAS_CASE_NAME} PROPERTIES TIMEOUT 2400)
endforeach()

# The Shape-DNA module against closed forms (src/ShapeDNA): exact reference
# integrals for P1-P3, the dense eigensolver against Eigen's, the Cholesky's
# solves, the rectangle and disk spectra (and the convergence orders 2, 4, 6),
# the square's exact double eigenvalues on a mesh with all its symmetries, the
# heat trace handing back area, perimeter and the corner term, invariance under
# similarity and remeshing, Lanczos against two dense solvers, and the
# spectrum bit-identical on one thread and on all of them.
add_test(NAME ShapeDNA_SelfTest COMMAND TestShapeDNA --selftest)
set_tests_properties(ShapeDNA_SelfTest PROPERTIES TIMEOUT 600)

# The default run (cubic elements, 50 eigenvalues, Dirichlet) on two corpus
# models: geom012, a plate with ten holes, and multimat/geom001, where the
# interface is ignored and the DNA is the unit square's. Exit code: convergence.
foreach(SHAPEDNA_CASE singlemat/geom012 multimat/geom001)
    string(REPLACE "/" "_" SHAPEDNA_CASE_NAME ${SHAPEDNA_CASE})
    add_test(NAME ShapeDNA_${SHAPEDNA_CASE_NAME}
             COMMAND TestShapeDNA ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/${SHAPEDNA_CASE}.obj --heat)
    set_tests_properties(ShapeDNA_${SHAPEDNA_CASE_NAME} PROPERTIES TIMEOUT 600)
endforeach()

# The TMOP smoother and the node-mobility layer under it, on meshes built in
# memory: a distorted grid whose boundary is straight (so a shape metric has to
# recover it exactly), an annulus whose two rims are closed feature loops with
# no corner on them (the case the closed-loop anchor exists for), a
# two-material grid (an interface is a feature curve that is not on dS), and a
# grid with one node dragged across its own one-ring (the untangler). It also
# differentiates the energy by hand against the analytic gradient and Hessian,
# and checks that one thread and many threads give bit-identical meshes -- the
# property the colouring is there to buy.
add_test(NAME TMOP_SelfTest COMMAND TestTMOP --selftest)

# TMOP end to end, on a mesh the pipeline actually produced rather than one
# built for it. geom012 is the model to run it on: a dozen voids, so 337 of its
# 1014 nodes are on a feature and slide rather than move freely, and it is the
# one case in the corpus where minimising the *average* metric (exponent 1)
# drives the worst scaled Jacobian down from 0.625 to 0.512 -- below where it
# started. At the default exponent of 2 it goes up, to 0.731, and the run
# asserts that the worst element is no worse than it was.
add_test(NAME MERIDIAN_TMOP_singlemat_geom012
         COMMAND TestMERIDIAN ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/singlemat/geom012.obj
                 --tmop 500)
set_tests_properties(MERIDIAN_TMOP_singlemat_geom012 PROPERTIES TIMEOUT 2400)

# Optional Qt6 OpenGL viewer (replaces the old GLFW-based viewer)
option(BUILD_OPENGL_VIEWER "Build the Qt6 OpenGL mesh/field viewer" OFF)
if (BUILD_OPENGL_VIEWER)
    find_package(OpenGL REQUIRED)
    # Svg is for the figure export: the GL feedback capture comes out as an SVG
    # document, and Qt paints that onto a PDF page.
    find_package(Qt6 REQUIRED COMPONENTS Core Widgets OpenGL OpenGLWidgets Svg)

    # Qt6 requires these global settings
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)

    add_executable(Viewer
        src/viewer/ViewerMain.cxx
        src/viewer/CrossGenWidget.cxx
        src/viewer/Geometry.cxx
        src/viewer/Render.cxx
        src/viewer/Interaction.cxx
        src/viewer/Export.cxx
    )
    target_include_directories(Viewer PRIVATE src)
    target_link_libraries(Viewer PRIVATE
        PolyVector
        OpenGL::GL
        Qt6::Core
        Qt6::Widgets
        Qt6::OpenGL
        Qt6::OpenGLWidgets
        Qt6::Svg
    )
    target_compile_definitions(Viewer PRIVATE GL_SILENCE_DEPRECATION)
    target_compile_features(Viewer PRIVATE cxx_std_17)
    message(STATUS "Building Qt6 OpenGL viewer")
endif()
add_executable(TopoDiag src/utils/TopoDiag.cxx)
target_link_libraries(TopoDiag PRIVATE PolyVector)

# The Python API (src/python): `crossgen`, a CPython extension module written
# against Python.h, over PolyVector.
#
#     PYTHONPATH=build/python python3 -c "import crossgen; help(crossgen)"
#
# It is built for the interpreter FindPython picks (the first python3 on PATH;
# -DPython3_EXECUTABLE=... chooses another) and imports only into that Python's
# minor version. Development.Module is the headers-only part of Python: an
# extension module takes its Python symbols from the interpreter that loads it
# (on macOS through -undefined dynamic_lookup, which Python3::Module carries),
# and linking one against libpython instead would give a process two
# interpreters' worth of state the moment it is loaded into a Python other than
# the one it was built against. numpy is not a build dependency: the module
# imports it at run time to hand back arrays.
option(BUILD_PYTHON_MODULE "Build the crossgen CPython extension module" ON)
if (BUILD_PYTHON_MODULE)
    if (CMAKE_VERSION VERSION_LESS 3.18)
        message(FATAL_ERROR "BUILD_PYTHON_MODULE needs CMake 3.18 or later "
                            "(FindPython's Development.Module); or pass -DBUILD_PYTHON_MODULE=OFF")
    endif()
    # The python3 on PATH -- the one that will `import crossgen` -- rather than
    # whichever framework FindPython reaches first on macOS, which on a machine
    # with Homebrew's Python beside python.org's is the wrong one.
    if (NOT DEFINED Python3_FIND_FRAMEWORK)
        set(Python3_FIND_FRAMEWORK LAST)
    endif()
    if (NOT DEFINED Python3_FIND_UNVERSIONED_NAMES)
        set(Python3_FIND_UNVERSIONED_NAMES FIRST)
    endif()
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
    message(STATUS "Python module for ${Python3_EXECUTABLE} (${Python3_VERSION})")

    Python3_add_library(crossgen MODULE WITH_SOABI
        src/python/PyOptions.hxx
        src/python/PyOptions.cxx
        src/python/PyUtil.hxx
        src/python/PyUtil.cxx
        src/python/PyTables.hxx
        src/python/PyTables.cxx
        src/python/Methods.hxx
        src/python/Methods.cxx
        src/python/MethodZIPLINE.cxx
        src/python/MethodUMBER.cxx
        src/python/MethodLayout.cxx
        src/python/MethodATLAS.cxx
        src/python/CrossGenModule.cxx
    )
    target_link_libraries(crossgen PRIVATE PolyVector)
    # Its own directory, so that PYTHONPATH can point at the module and at
    # nothing else of the build tree.
    set_target_properties(crossgen PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python)

    # A build cannot change the PYTHONPATH of the shell that started it, so this
    # is the other way to make `import crossgen` work without one: a one-line
    # .pth file in the site-packages of the interpreter the module was built
    # for, naming build/python. Python reads it at start-up in every shell,
    # IDE and notebook, and it follows rebuilds because it names the directory
    # and not the .so. The site directory is the user's when the interpreter
    # has one, and the environment's own inside a venv, where the user site is
    # switched off. It is written at configure time and only when asked for,
    # since it puts a file outside the repository.
    option(CROSSGEN_PYTHON_DEV_PTH "Make the built crossgen module importable from this Python (writes crossgen-dev.pth into its site-packages)" OFF)
    if (CROSSGEN_PYTHON_DEV_PTH)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c
                "import site, sysconfig; print(site.getusersitepackages() if site.ENABLE_USER_SITE else sysconfig.get_paths()['purelib'])"
            OUTPUT_VARIABLE CROSSGEN_PTH_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE CROSSGEN_PTH_RC)
        if (NOT CROSSGEN_PTH_RC EQUAL 0 OR CROSSGEN_PTH_DIR STREQUAL "")
            message(FATAL_ERROR "CROSSGEN_PYTHON_DEV_PTH: could not ask ${Python3_EXECUTABLE} for its site-packages")
        endif()
        file(MAKE_DIRECTORY "${CROSSGEN_PTH_DIR}")
        file(WRITE "${CROSSGEN_PTH_DIR}/crossgen-dev.pth" "${CMAKE_BINARY_DIR}/python\n")
        message(STATUS "crossgen importable: ${CROSSGEN_PTH_DIR}/crossgen-dev.pth -> ${CMAKE_BINARY_DIR}/python "
                       "(delete that file to undo)")
    endif()

    # Every Mesh method end to end on one model -- a block decomposition, a
    # mesh on it and a few smoothing sweeps -- plus the Shape-DNA, the array
    # round trip, and the keyword errors. See the script's header.
    add_test(NAME Python_Module
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/python/test_crossgen.py
                     ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/singlemat/geom001.obj)
    set_tests_properties(Python_Module PROPERTIES
                         ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python"
                         TIMEOUT 600)
endif()

# The Tier-1 experiments of the paper (Sec. 6): E1 to E5.
add_subdirectory(paper_tests)
