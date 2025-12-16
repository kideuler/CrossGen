# Try to locate Gmsh library and headers
# Provides: Gmsh::gmsh target

find_path(GMSH_INCLUDE_DIR gmsh.h
  HINTS $ENV{GMSH_ROOT}/include /usr/include /usr/local/include /opt/homebrew/include
)
find_library(GMSH_LIBRARY NAMES gmsh
  HINTS $ENV{GMSH_ROOT}/lib /usr/lib /usr/local/lib /opt/homebrew/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh REQUIRED_VARS GMSH_INCLUDE_DIR GMSH_LIBRARY)

if(GMSH_FOUND)
  add_library(Gmsh::gmsh UNKNOWN IMPORTED)
  set_target_properties(Gmsh::gmsh PROPERTIES
    IMPORTED_LOCATION "${GMSH_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
  )
endif()
