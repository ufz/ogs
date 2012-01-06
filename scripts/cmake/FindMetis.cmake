#
# Find the METIS includes and libraries
#
# METIS is an library that implements a variety of algorithms for
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of
# sparse matrices. It can be found at:
# http://glaros.dtc.umn.edu/gkhome/metis/metis/download
#
# METIS_INCLUDE_DIR - where to find autopack.h
# METIS_LIBRARIES   - List of fully qualified libraries to link against.
# METIS_FOUND       - Do not attempt to use if "no" or undefined.

FIND_PATH(METIS_INCLUDE_DIR metis.h
	/usr/include/metis
	$ENV{HOME}/include/
)

FIND_LIBRARY(METIS_LIBRARY metis
	/usr/lib
	$ENV{HOME}/lib/
)

IF(METIS_INCLUDE_DIR)
  IF(METIS_LIBRARY)
    SET( METIS_LIBRARIES ${METIS_LIBRARY})
    SET( METIS_FOUND "YES" )
  ENDIF(METIS_LIBRARY)
ENDIF(METIS_INCLUDE_DIR)
