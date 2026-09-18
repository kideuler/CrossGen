# OpenCASCADE is src/geom's business and nobody else's.
#
#   cmake -DROOT=<repository> -DOCC_INCLUDE=<OpenCASCADE include dir> -P CheckOpenCascadeConfinement.cmake
#
# Fails, naming the file and the header, if any source outside src/geom includes
# a header that exists in OpenCASCADE's include directory, or includes the
# bridge header geom/detail/Occ.hxx. The build already refuses such a file --
# the kernel's include directories are private to the CrossGenGeom target --
# but only once it is compiled, and only on a machine where the kernel is not
# also on a default include path. This says it on the sources, everywhere.

if (NOT ROOT OR NOT OCC_INCLUDE)
    message(FATAL_ERROR "usage: cmake -DROOT=<repo> -DOCC_INCLUDE=<dir> -P ${CMAKE_CURRENT_LIST_FILE}")
endif()
if (NOT EXISTS "${OCC_INCLUDE}/Standard.hxx")
    message(FATAL_ERROR "${OCC_INCLUDE} does not look like OpenCASCADE's include directory")
endif()

file(GLOB_RECURSE sources
     "${ROOT}/src/*.cxx" "${ROOT}/src/*.hxx" "${ROOT}/src/*.cpp" "${ROOT}/src/*.hpp"
     "${ROOT}/src/*.h" "${ROOT}/paper_tests/*.cxx" "${ROOT}/paper_tests/*.hxx")

set(geom_dir "${ROOT}/src/geom/")
string(LENGTH "${geom_dir}" geom_len)
set(offenders "")
set(checked 0)
foreach(source ${sources})
    string(SUBSTRING "${source}" 0 ${geom_len} prefix)
    if (prefix STREQUAL geom_dir)
        continue()
    endif()
    math(EXPR checked "${checked} + 1")
    file(STRINGS "${source}" lines REGEX "^[ \t]*#[ \t]*include[ \t]*[<\"][^>\"]+[>\"]")
    foreach(line ${lines})
        string(REGEX REPLACE "^[ \t]*#[ \t]*include[ \t]*[<\"]([^>\"]+)[>\"].*$" "\\1" header "${line}")
        if (header MATCHES "^geom/detail/")
            list(APPEND offenders "${source}: ${header} (the bridge to the kernel is internal to src/geom)")
        elseif (NOT header MATCHES "/" AND EXISTS "${OCC_INCLUDE}/${header}")
            list(APPEND offenders "${source}: ${header}")
        endif()
    endforeach()
endforeach()

if (offenders)
    string(REPLACE ";" "\n  " report "${offenders}")
    message(FATAL_ERROR "OpenCASCADE is used outside src/geom -- add what is needed to src/geom instead:\n  ${report}")
endif()
message(STATUS "${checked} source file(s) outside src/geom include nothing of OpenCASCADE")
