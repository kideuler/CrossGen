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
set(COMISO_ENABLE_DEFAULT OFF CACHE BOOL "" FORCE)
add_subdirectory(extern/CoMiSo)
message(STATUS "Using local CoMiSo from: ${CMAKE_CURRENT_SOURCE_DIR}/extern/CoMiSo")

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
    src/polyvector/PolyVectors.hxx
    src/polyvector/PolyVectors.cxx
    src/SIPG/SIPG.hxx
    src/SIPG/SIPG.cxx
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
target_link_libraries(PolyVector PUBLIC CoMISo::CoMISo)

# Small utility to load and print mesh data
add_executable(LoadMesh src/utils/LoadMesh.cxx)
target_link_libraries(LoadMesh PRIVATE PolyVector)

# Utility to trace separatrices on a mesh
add_executable(TraceMesh src/utils/TraceMesh.cxx)
target_link_libraries(TraceMesh PRIVATE PolyVector)

# Utility to cut a mesh into a disk and run sanity checks
add_executable(TestCutMesh src/utils/TestCutMesh.cxx)
target_link_libraries(TestCutMesh PRIVATE PolyVector)

# Test MIQ parametrization quality
add_executable(TestMIQ src/utils/TestMIQ.cxx)
target_link_libraries(TestMIQ PRIVATE PolyVector)

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

# Test p=0 DG/SIPG MBO cross-field solver
add_executable(TestSIPG src/utils/TestSIPG.cxx)
target_link_libraries(TestSIPG PRIVATE PolyVector)

# Test the UMBER frame field optimization (Wang et al. 2022, Secs. 4.1-4.2)
add_executable(TestUMBER src/utils/TestUMBER.cxx)
target_link_libraries(TestUMBER PRIVATE PolyVector)

# Test MERIDIAN Stages 1-6 (Shepherd, Gu and Hughes 2022, Secs. 3.1-3.3):
# cone singularities, cutting graph, Gauss-Bonnet, discrete Ricci flow
add_executable(TestMERIDIAN src/utils/TestMERIDIAN.cxx)
target_link_libraries(TestMERIDIAN PRIVATE PolyVector)

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

# Custom target to run TestMIQ on all meshes
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

# Custom target to run TestSIPG on all meshes
set(SIPG_COMMANDS "")
foreach(MESH_FILE ${MESH_OBJ_FILES})
    get_filename_component(MESH_NAME ${MESH_FILE} NAME)
    list(APPEND SIPG_COMMANDS 
        COMMAND ${CMAKE_COMMAND} -E echo "Processing ${MESH_NAME}..."
        COMMAND $<TARGET_FILE:TestSIPG> ${MESH_FILE} || ${CMAKE_COMMAND} -E true
    )
endforeach()

add_custom_target(sipg
    ${SIPG_COMMANDS}
    DEPENDS TestSIPG
    COMMENT "Running TestSIPG on all meshes in data/meshes"
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

# End-to-end: Stages 1 to 6 on a handful of models that between them exercise
# every part of the later stages -- a plain disk with one interior cone
# (geom001), a model whose cutting graph carries several different holonomies
# (geom003), one with voids and a dozen connectivity constraints (geom012), and
# a two-material model, where the interface between the physical surfaces is a
# feature chain and E3 has something to hold (multimat geom001). Each asserts
# Q1 to Q5 on the map that comes out; the whole corpus is `make meridian`.
foreach(MERIDIAN_CASE singlemat/geom001 singlemat/geom003 singlemat/geom012 multimat/geom001)
    string(REPLACE "/" "_" MERIDIAN_CASE_NAME ${MERIDIAN_CASE})
    add_test(NAME MERIDIAN_Layout_${MERIDIAN_CASE_NAME}
             COMMAND TestMERIDIAN ${CMAKE_CURRENT_SOURCE_DIR}/data/meshes/${MERIDIAN_CASE}.obj)
endforeach()
# Optional Qt6 OpenGL viewer (replaces the old GLFW-based viewer)
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