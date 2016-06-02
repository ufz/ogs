MESSAGE( STATUS "The METIS package is copyrighted by the Regents of the University of Minnesota." )
MESSAGE( STATUS "Please read the license of the METIS package carefully before you use the METIS." )

add_definitions(-DUSE_GKREGEX)

set(GKLIB_PATH "${METIS_PATH}/GKlib" CACHE PATH "path to GKlib")
set(SHARED FALSE CACHE BOOL "build a shared library")

if(MSVC)
	set(METIS_INSTALL FALSE)
else()
	set(METIS_INSTALL TRUE)
endif()

# Configure libmetis library.
if(SHARED)
	set(METIS_LIBRARY_TYPE SHARED)
else()
	set(METIS_LIBRARY_TYPE STATIC)
endif(SHARED)

include(${GKLIB_PATH}/GKlibSystem.cmake)
# Add include directories.
include_directories(${GKLIB_PATH})
include_directories(${METIS_PATH}/include)
# Recursively look for CMakeLists.txt in subdirs.
add_subdirectory("${METIS_PATH}/include")
add_subdirectory("${METIS_PATH}/libmetis")
