# - Try to find SCIP (https://scipopt.org)
#
# Once done, this will define:
#   SCIP_FOUND
#   SCIP_INCLUDE_DIRS
#   SCIP_LIBRARIES
#   SCIP::SCIP (imported target)
#
# Hints:
#   -DSCIP_ROOT=<prefix>
#   -DSCIP_DIR=<prefix>

include(FindPackageHandleStandardArgs)

set(_SCIP_HINT_DIRS "")
if(DEFINED SCIP_ROOT)
  list(APPEND _SCIP_HINT_DIRS "${SCIP_ROOT}")
endif()
if(DEFINED SCIP_DIR)
  list(APPEND _SCIP_HINT_DIRS "${SCIP_DIR}")
endif()
if(DEFINED ENV{SCIP_ROOT})
  list(APPEND _SCIP_HINT_DIRS "$ENV{SCIP_ROOT}")
endif()
if(DEFINED ENV{SCIP_DIR})
  list(APPEND _SCIP_HINT_DIRS "$ENV{SCIP_DIR}")
endif()

find_path(SCIP_INCLUDE_DIR
  NAMES scip/scip.h
  HINTS ${_SCIP_HINT_DIRS}
  PATH_SUFFIXES include
)

find_library(SCIP_LIBRARY
  NAMES scip libscip
  HINTS ${_SCIP_HINT_DIRS}
  PATH_SUFFIXES lib lib64
)

set(SCIP_INCLUDE_DIRS ${SCIP_INCLUDE_DIR})
set(SCIP_LIBRARIES ${SCIP_LIBRARY})

find_package_handle_standard_args(SCIP
  REQUIRED_VARS SCIP_INCLUDE_DIR SCIP_LIBRARY
)

if(SCIP_FOUND AND NOT TARGET SCIP::SCIP)
  add_library(SCIP::SCIP UNKNOWN IMPORTED)
  set_target_properties(SCIP::SCIP PROPERTIES
    IMPORTED_LOCATION "${SCIP_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${SCIP_INCLUDE_DIR}"
  )
endif()
