message(
    STATUS
        "The METIS package is copyrighted by the Regents of the University of Minnesota."
        "Please read the license of the METIS package carefully before you use METIS."
)

# cmake-lint: disable=C0103

# Metis library
file(GLOB _metis_sources ${metis_SOURCE_DIR}/libmetis/*.c)
file(GLOB GKlib_sources ${GKlib_SOURCE_DIR}/*.c)
if(WIN32)
    set(_metis_static STATIC)
endif()
ogs_add_library(ogs_metis ${_metis_static} ${GKlib_sources} ${_metis_sources})
target_compile_definitions(
    ogs_metis PUBLIC USE_GKREGEX
                     $<$<CXX_COMPILER_ID:MSVC>:__thread=__declspec\(thread\)>
)
target_include_directories(
    ogs_metis PUBLIC ${GKlib_SOURCE_DIR} ${metis_SOURCE_DIR}/include
                     ${metis_SOURCE_DIR}/libmetis
)
if(OpenMP_FOUND)
    target_link_libraries(ogs_metis PUBLIC OpenMP::OpenMP_C)
endif()
if(UNIX)
    target_link_libraries(ogs_metis PUBLIC m)
endif()
target_compile_options(
    ogs_metis PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
                      $<$<CXX_COMPILER_ID:MSVC>:/W0>
)
target_compile_definitions(ogs_metis PUBLIC IDXTYPEWIDTH=64 REALTYPEWIDTH=32)
set_property(TARGET ogs_metis PROPERTY C_STANDARD 90)

# mpmetis binary
set(_mpmetis_sources
    ${metis_SOURCE_DIR}/programs/mpmetis.c
    ${metis_SOURCE_DIR}/programs/cmdline_mpmetis.c
    ${metis_SOURCE_DIR}/programs/io.c ${metis_SOURCE_DIR}/programs/stat.c
)
ogs_add_executable(mpmetis ${_mpmetis_sources})
target_compile_definitions(mpmetis PRIVATE SVNINFO="")
target_link_libraries(mpmetis ogs_metis)
target_compile_options(
    mpmetis PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
                    $<$<CXX_COMPILER_ID:MSVC>:/W0>
)
set_property(TARGET mpmetis PROPERTY C_STANDARD 90)

install(TARGETS mpmetis RUNTIME DESTINATION bin)
