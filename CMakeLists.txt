cmake_minimum_required(VERSION 3.15)
project(CrossGen LANGUAGES CXX)

# Build options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Eigen (header-only)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "Eigen3 include dir: ${EIGEN3_INCLUDE_DIR}")

option(BUILD_MESH2GMSH "Build the Mesh2Dgmsh utility" ON)

if (BUILD_MESH2GMSH)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
    find_package(Gmsh REQUIRED)
    message(STATUS "Gmsh found: ${GMSH_LIBRARY}")
    add_executable(Mesh2Dgmsh src/Mesh2Dgmsh.cxx)
    target_link_libraries(Mesh2Dgmsh PRIVATE Gmsh::gmsh)
else()
    message(STATUS "Skipping Mesh2Dgmsh utility")
endif()

add_library(PolyVector SHARED
    src/polyvec/VertexTriangleCSR.hxx
    src/polyvec/VertexTriangleCSR.cxx
    src/polyvec/Mesh.cxx
    src/polyvec/Mesh.hxx
    src/polyvec/PolyVectors.hxx
    src/polyvec/PolyVectors.cxx
)
target_include_directories(PolyVector PUBLIC src/polyvec)
target_link_libraries(PolyVector PUBLIC Eigen3::Eigen)

# Small utility to load and print mesh data
add_executable(LoadMesh src/LoadMesh.cxx)
target_link_libraries(LoadMesh PRIVATE PolyVector)

# Optional OpenGL viewer
option(BUILD_OPENGL_VIEWER "Build the OpenGL mesh/field viewer" OFF)
if (BUILD_OPENGL_VIEWER)
    if(APPLE)
        find_package(OpenGL REQUIRED)
        # Try multiple glfw package names/targets for robustness on macOS/Homebrew
        set(GLFW_TARGET "")
        # First try config package (defines target `glfw` or `glfw::glfw`)
        find_package(glfw3 CONFIG QUIET)
        if (TARGET glfw)
            set(GLFW_TARGET glfw)
        elseif (TARGET glfw::glfw)
            set(GLFW_TARGET glfw::glfw)
        endif()

        # Fallback to module mode (may also define `glfw` or `glfw3::glfw`)
        if (NOT GLFW_TARGET)
            find_package(glfw3 QUIET)
            if (TARGET glfw)
                set(GLFW_TARGET glfw)
            elseif (TARGET glfw3::glfw)
                set(GLFW_TARGET glfw3::glfw)
            endif()
        endif()

        if (GLFW_TARGET)
            message(STATUS "Building OpenGL viewer with target ${GLFW_TARGET}")
            add_executable(Viewer src/Viewer.cxx)
            target_link_libraries(Viewer PRIVATE PolyVector OpenGL::GL ${GLFW_TARGET})
            target_compile_definitions(Viewer PRIVATE GL_SILENCE_DEPRECATION GLFW_INCLUDE_NONE)
            # Prefer a modern C++ standard for the viewer too
            target_compile_features(Viewer PRIVATE cxx_std_17)
        else()
            message(WARNING "glfw3 not found. The OpenGL viewer will not be built. Install via Homebrew: brew install glfw")
        endif()
    else()
        message(WARNING "BUILD_OPENGL_VIEWER is currently supported only on macOS in this project setup.")
    endif()
endif()