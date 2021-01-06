# - Try to find libgeotiff
#
# Once done, this will define
#
#  GEOTIFF_FOUND
#  GEOTIFF_INCLUDE_DIRS
#  GEOTIFF_LIBRARIES

###
# Dependencies
###
set(_deps_libs)
set(_deps_includes)
set(_deps_check)

find_path( libgeotiff_INCLUDE_DIR geotiff.h)
find_library(libgeotiff_LIBRARY geotiff)

find_path( xtiff_INCLUDE_DIR xtiffio.h)
if(MSVC)
    find_library(xtiff_LIBRARY xtiff)
    list(APPEND _deps_libs ${xtiff_LIBRARY})
endif()

find_package(TIFF)

list(APPEND _deps_libs ${TIFF_LIBRARIES})
list(APPEND _deps_includes ${TIFF_INCLUDE_DIRS})
list(APPEND _deps_check TIFF_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GEOTIFF
    REQUIRED_VARS
    libgeotiff_LIBRARY
    libgeotiff_INCLUDE_DIR
    xtiff_INCLUDE_DIR
    ${_deps_check}
)

if(GEOTIFF_FOUND)
    set(GEOTIFF_INCLUDE_DIRS ${libgeotiff_INCLUDE_DIR} ${xtiff_INCLUDE_DIR} ${_deps_includes})
    set(GEOTIFF_LIBRARIES ${libgeotiff_LIBRARY} ${_deps_libs})
endif()
