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

# CoMiSo library - use local submodule
#
# CoMiSo is only needed for MIQ (mixed-integer quadrangulation); it is a large
# dependency to build and most work in this repo never touches it. Turn this
# off to skip building it entirely; code that depends on it is compiled out
# via the CROSSGEN_WITH_COMISO define (see MIQ.hxx / MIQ.cxx).
option(CROSSGEN_ENABLE_COMISO "Build the bundled CoMiSo solver and enable MIQ" OFF)

if (CROSSGEN_ENABLE_COMISO)
    set(COMISO_ENABLE_DEFAULT OFF CACHE BOOL "" FORCE)
    add_subdirectory(extern/CoMiSo)
    message(STATUS "Using local CoMiSo from: ${CMAKE_CURRENT_SOURCE_DIR}/extern/CoMiSo")
else()
    message(STATUS "CROSSGEN_ENABLE_COMISO is OFF: skipping CoMiSo build (MIQ disabled)")
endif()

# OpenMP
#
# Used by the TMOP smoother (src/mesh/TMOP.cxx) and by nothing else so far. The
# smoother sweeps in colour order, so its answer does not depend on the thread
# count and a build without OpenMP is a correct build, only a slower one --
# every pragma is compiled out through _OPENMP.
#
# AppleClang is the awkward case: it understands the pragmas but ships no
# runtime, so plain find_package(OpenMP) fails on a Mac that has Homebrew's
# libomp installed and working. Fall back to wiring libomp up by hand there.
option(CROSSGEN_ENABLE_OPENMP "Use OpenMP to parallelise the TMOP smoother" ON)

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
    message(STATUS "No OpenMP: TMOP will run single-threaded")
endif()

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
    src/polyvector/PolyVectors.hxx
    src/polyvector/PolyVectors.cxx
    src/dualmbo/DualMBO.hxx
    src/dualmbo/DualMBO.cxx
    src/Parameterization/CutMesh.hxx
    src/Parameterization/CutMesh.cxx
    src/Parameterization/HarmonicCut.hxx
    src/Parameterization/HarmonicCut.cxx
    src/Parameterization/MIQ.hxx
    src/Parameterization/MIQ.cxx
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
    src/UMBER/Polysquare.hxx
    src/UMBER/Polysquare.cxx
    src/UMBER/MotorcycleGraph.hxx
    src/UMBER/MotorcycleGraph.cxx
    src/UMBER/BlockLayout.hxx
    src/UMBER/BlockLayout.cxx
    src/UMBER/ChordCollapse.hxx
    src/UMBER/ChordCollapse.cxx
    src/tracing/FieldTracer.hxx
    src/tracing/FieldTracer.cxx
    src/tracing/SeparatrixTrace.hxx
    src/tracing/SeparatrixTrace.cxx
    src/tracing/QuadLayout.hxx
    src/tracing/QuadLayout.cxx
    src/tracing/PartitionSimplify.hxx
    src/tracing/PartitionSimplify.cxx
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
target_link_libraries(PolyVector PUBLIC Eigen3::Eigen)
if (CROSSGEN_OPENMP_TARGET)
    target_link_libraries(PolyVector PUBLIC ${CROSSGEN_OPENMP_TARGET})
endif()
if (CROSSGEN_ENABLE_COMISO)
    target_link_libraries(PolyVector PUBLIC CoMISo::CoMISo)
    target_compile_definitions(PolyVector PUBLIC CROSSGEN_WITH_COMISO)
endif()

# Small utility to load and print mesh data
add_executable(LoadMesh src/utils/LoadMesh.cxx)
target_link_libraries(LoadMesh PRIVATE PolyVector)

# Utility to trace separatrices on a mesh
add_executable(TraceMesh src/utils/TraceMesh.cxx)
target_link_libraries(TraceMesh PRIVATE PolyVector)

# Utility to cut a mesh into a disk and run sanity checks
add_executable(TestCutMesh src/utils/TestCutMesh.cxx)
target_link_libraries(TestCutMesh PRIVATE PolyVector)

# Test MIQ parametrization quality (needs CoMiSo)
if (CROSSGEN_ENABLE_COMISO)
    add_executable(TestMIQ src/utils/TestMIQ.cxx)
    target_link_libraries(TestMIQ PRIVATE PolyVector)
endif()

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

# Test the UMBER frame field optimization (Wang et al. 2022, Secs. 4.1-4.2)
add_executable(TestUMBER src/utils/TestUMBER.cxx)
target_link_libraries(TestUMBER PRIVATE PolyVector)

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

# Custom target to run TestMIQ on all meshes (needs CoMiSo)
if (CROSSGEN_ENABLE_COMISO)
    set(MIQ_COMMANDS "")
    foreach(MESH_FILE ${MESH_OBJ_FILES})
        get_filename_component(MESH_NAME ${MESH_FILE} NAME)
        list(APPEND MIQ_COMMANDS
            COMMAND ${CMAKE_COMMAND} -E echo "Processing ${MESH_NAME}..."
            COMMAND $<TARGET_FILE:TestMIQ> ${MESH_FILE} || ${CMAKE_COMMAND} -E true
        )
    endforeach()

    add_custom_target(miq
        ${MIQ_COMMANDS}
        DEPENDS TestMIQ
        COMMENT "Running TestMIQ on all meshes in data/meshes"
        VERBATIM
    )
endif()

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
# graph whose faces tile the model exactly with nothing left dangling. How many
# of the components come out four-sided is printed but not asserted -- that is
# the part Sec. 4 of the paper is for.
add_test(NAME Tracing_LayoutsOnAllMeshes
         COMMAND TraceMesh --no-vtu ${MESH_OBJ_FILES})

# Discrete Ricci flow (Shepherd et al. Sec. 3.2.1) against prescribed cone sets
# with far more curvature to move than anything in data/meshes asks for, which
# is what puts the weighted-Delaunay flipping under load. Checks the cone angle
# sums directly and the Hessian against a finite difference.
add_test(NAME MERIDIAN_RicciFlowSelfTest COMMAND TestMERIDIAN --selftest)

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
#
# Builds without CoMiSo too: the UV-mesh split-screen view (phase 5, MIQ) is
# compiled out via CROSSGEN_WITH_COMISO and that phase is simply unreachable.
option(BUILD_OPENGL_VIEWER "Build the Qt6 OpenGL mesh/field viewer" OFF)
if (BUILD_OPENGL_VIEWER)
    find_package(OpenGL REQUIRED)
    find_package(Qt6 REQUIRED COMPONENTS Core Widgets OpenGL OpenGLWidgets)

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
    )
    target_compile_definitions(Viewer PRIVATE GL_SILENCE_DEPRECATION)
    target_compile_features(Viewer PRIVATE cxx_std_17)
    message(STATUS "Building Qt6 OpenGL viewer")
endif()
add_executable(TopoDiag src/utils/TopoDiag.cxx)
target_link_libraries(TopoDiag PRIVATE PolyVector)

# The Tier-1 experiments of the paper (Sec. 6): E1 to E5.
add_subdirectory(paper_tests)
