#
# Find the METIS includes and libraries
#
#    METIS_INCLUDE_DIR - where to find metis.h
#    METIS_LIBRARIES   - List of fully qualified libraries to link against.
#    METIS_FOUND       - Do not attempt to use if "no" or undefined.
#
# METIS is an library that implements a variety of algorithms for
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of
# sparse matrices. It can be found at:
# http://glaros.dtc.umn.edu/gkhome/metis/metis/download
#

FIND_PATH(METIS_INCLUDE_DIR metis.h
	/usr/include/metis
	$ENV{HOME}/include/
	${CMAKE_SOURCE_DIR}/../libs/include
	${OGS_LIBS_DIR_FOUND}/include
)

FIND_LIBRARY(METIS_LIBRARY metis
	/usr/lib
	$ENV{HOME}/lib/
	${CMAKE_SOURCE_DIR}/../libs/lib
	${OGS_LIBS_DIR_FOUND}/lib
)

SET(METIS_LIBRARIES ${METIS_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Metis DEFAULT_MSG METIS_LIBRARY METIS_INCLUDE_DIR)
