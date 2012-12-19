# - Try to find LIS
# Once done, this will define
#
#  LIS_FOUND
#  LIS_INCLUDE_DIRS
#  LIS_LIBRARIES

set(LIS_ROOT_DIR
    "${LIS_ROOT_DIR}"
    CACHE
    PATH
    "Directory to search for Lis library")
	
find_path( LIS_INCLUDE_DIR
	NAMES lis.h
	HINTS
	${LIS_ROOT_DIR}/include
    /usr/include/lis
    $ENV{HOME}/include/
	)

find_library(LIS_LIBRARY
    NAMES lis
    HINTS 
    ${LIS_ROOT_DIR}/lib
    /usr/lib
    $ENV{HOME}/lib/
    ) 

SET(LIS_LIBRARIES ${LIS_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIS DEFAULT_MSG LIS_LIBRARY LIS_INCLUDE_DIR)

