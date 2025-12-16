cmake_minimum_required(VERSION 3.15)
project(CrossGen LANGUAGES CXX)

# Build options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Gmsh REQUIRED)
message(STATUS "Gmsh found: ${GMSH_LIBRARY}")


add_executable(Mesh2Dgmsh src/Mesh2Dgmsh.cpp)
target_link_libraries(Mesh2Dgmsh PRIVATE Gmsh::gmsh)