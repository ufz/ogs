message( STATUS "The METIS package is copyrighted by the Regents of the University of Minnesota." )
message( STATUS "Please read the license of the METIS package carefully before you use the METIS." )

set(METIS_PATH ${PROJECT_SOURCE_DIR}/ThirdParty/metis)
add_definitions(-DUSE_GKREGEX)
set(GKLIB_PATH "${METIS_PATH}/GKlib" CACHE PATH "path to GKlib")

if(SHARED)
	set(METIS_LIBRARY_TYPE SHARED)
else()
	set(METIS_LIBRARY_TYPE STATIC)
endif(SHARED)

include(${GKLIB_PATH}/GKlibSystem.cmake)
include_directories(${GKLIB_PATH})
include_directories(${METIS_PATH}/include)

# From ${METIS_PATH}/libmetis/CMakeLists.txt
# Removed linking to conan
# Add this directory for internal users.
include_directories(BEFORE ${METIS_PATH}/libmetis)
# Find sources.
file(GLOB metis_sources ${METIS_PATH}/libmetis/*.c)
# Build libmetis.
add_library(metis ${GKlib_sources} ${metis_sources})
if(BUILD_SHARED_LIBS)
    install(TARGETS metis LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(UNIX)
  target_link_libraries(metis m)
elseif(MSVC)
  include(GenerateExportHeader)
  generate_export_header(metis)
endif()

if(BUILD_SHARED_LIBS)
    install(TARGETS metis
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        COMPONENT ogs_partmesh
    )
endif()

## Compile mpmetis
add_definitions(-DSVNINFO="")
include_directories(${METIS_PATH}/libmetis)
include_directories(${METIS_PATH}/programs)
set(METIS_SOURCES
   ${METIS_PATH}/programs/mpmetis.c
   ${METIS_PATH}/programs/cmdline_mpmetis.c
   ${METIS_PATH}/programs/io.c
   ${METIS_PATH}/programs/stat.c
   )
add_executable(mpmetis ${METIS_SOURCES})
target_link_libraries(mpmetis metis)
install(TARGETS mpmetis RUNTIME DESTINATION bin COMPONENT ogs_partmesh)

# Disable warnings
if(MSVC)
    set_target_properties(metis mpmetis PROPERTIES COMPILE_FLAGS /W0)
else()
    set_target_properties(metis mpmetis PROPERTIES COMPILE_FLAGS -w)
endif()
