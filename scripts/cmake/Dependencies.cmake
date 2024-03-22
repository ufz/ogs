message(STATUS "┌─ Dependencies.cmake")
list(APPEND CMAKE_MESSAGE_INDENT "│    ")
set(CMAKE_FOLDER ThirdParty)

if(OGS_BUILD_TESTING)
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        VERSION ${ogs.minimum_version.gtest}
        GIT_TAG v${ogs.tested_version.gtest}
        OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt ON"
                "BUILD_SHARED_LIBS OFF"
        EXCLUDE_FROM_ALL YES
    )
    if(googletest_ADDED AND WIN32)
        target_compile_options(gtest PRIVATE /EHsc)
        target_compile_options(gmock PRIVATE /EHsc)
    endif()

    CPMAddPackage(
        NAME autocheck
        GITHUB_REPOSITORY ufz/autocheck
        GIT_TAG e388ecbb31c49fc2724c8d0436da313b6edca7fd
        DOWNLOAD_ONLY YES
    )
    if(autocheck_ADDED)
        add_library(autocheck INTERFACE IMPORTED)
        target_include_directories(
            autocheck SYSTEM INTERFACE ${autocheck_SOURCE_DIR}/include
        )
    endif()
endif()

if(GUIX_BUILD)
    find_package(spdlog REQUIRED)
else()
    CPMAddPackage(
        NAME spdlog
        GITHUB_REPOSITORY gabime/spdlog
        VERSION 1.12.0
        OPTIONS "BUILD_SHARED_LIBS OFF"
    )
endif()

if(GUIX_BUILD)
    add_library(tclap INTERFACE IMPORTED) # header-only, nothing else to do
else()
    CPMFindPackage(
        NAME tclap
        GITHUB_REPOSITORY ufz/tclap
        VERSION 1.2.4
        GIT_TAG 098dd0fe07a31618f3c2a9f8727bb01c8c5d61e2
        DOWNLOAD_ONLY YES
    )
    if(tclap_ADDED)
        add_library(tclap INTERFACE IMPORTED)
        target_include_directories(
            tclap SYSTEM INTERFACE ${tclap_SOURCE_DIR}/include
        )
    endif()
endif()

if(GUIX_BUILD)
    find_library(tet tet REQUIRED)
    find_program(TETGEN_EXECUTABLE tetgen REQUIRED)
    install(PROGRAMS ${TETGEN_EXECUTABLE} DESTINATION bin)
else()
    CPMAddPackage(
        NAME tetgen GITHUB_REPOSITORY ufz/tetgen
        GIT_TAG 213548f5bca1ec00269603703f0fec1272181587
    )
    if(tetgen_ADDED)
        install(PROGRAMS $<TARGET_FILE:tetgen> DESTINATION bin)
        list(APPEND DISABLE_WARNINGS_TARGETS tet tetgen)
    endif()
endif()

CPMFindPackage(
    NAME pybind11 GITHUB_REPOSITORY pybind/pybind11
    VERSION ${ogs.minimum_version.pybind11}
)

if(_build_chemistry_lib)
    if(GUIX_BUILD)
        find_library(
            IPhreeqc_LIBRARY NAMES iphreeqc IPhreeqc IPhreeqcrwdi REQUIRED
        )
        add_library(IPhreeqc INTERFACE IMPORTED)
        target_link_libraries(IPhreeqc INTERFACE ${IPhreeqc_LIBRARY})
    else()
        CPMAddPackage(
            NAME iphreeqc GITHUB_REPOSITORY ufz/iphreeqc GIT_TAG 3.5.0-1
        )
        if(iphreeqc_ADDED)
            target_include_directories(
                IPhreeqc
                PUBLIC ${iphreeqc_SOURCE_DIR}/src
                INTERFACE ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/common
                          ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/PhreeqcKeywords
                          ${iphreeqc_SOURCE_DIR}/src/phreeqcpp
            )
            list(APPEND DISABLE_WARNINGS_TARGETS IPhreeqc)
            if(BUILD_SHARED_LIBS)
                install(TARGETS IPhreeqc
                        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                )
            endif()
        endif()
    endif()
endif()

set(_eigen_version ${ogs.minimum_version.eigen})
set(_eigen_url
    https://gitlab.com/libeigen/eigen/-/archive/${_eigen_version}/eigen-${_eigen_version}.tar.gz
)
if(OGS_USE_EIGEN_UNSUPPORTED)
    set(_eigen_version 3.4.90)
    set(_eigen_url
        https://gitlab.com/libeigen/eigen/-/archive/${ogs.minimum_version.eigen-unsupported}/eigen-${ogs.minimum_version.eigen-unsupported}.tar.gz
    )
endif()

CPMFindPackage(
    NAME Eigen3
    # Error as in
    # https://gitlab.com/gitlab-com/gl-infra/reliability/-/issues/8475
    # GITLAB_REPOSITORY libeigen/eigen
    URL ${_eigen_url}
    VERSION ${_eigen_version}
    DOWNLOAD_ONLY YES
)
if(Eigen3_ADDED)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(
        Eigen3::Eigen SYSTEM INTERFACE ${Eigen3_SOURCE_DIR}
                                       ${OpenMP_CXX_INCLUDE_DIRS}
    )
endif()

if(OGS_USE_MFRONT)
    if(GUIX_BUILD)
        find_package(MFrontGenericInterface REQUIRED)
    else()
        set(CMAKE_REQUIRE_FIND_PACKAGE_TFEL TRUE)
        CPMAddPackage(
            NAME MGIS
            GITHUB_REPOSITORY thelfer/MFrontGenericInterfaceSupport
            GIT_TAG rliv-2.0
            OPTIONS "enable-doxygen-doc OFF" "enable-fortran-bindings OFF"
                    "enable-website OFF"
            EXCLUDE_FROM_ALL YES
        )
        if(MGIS_ADDED)
            list(APPEND DISABLE_WARNINGS_TARGETS MFrontGenericInterface)
        endif()
    endif()
endif()

CPMFindPackage(
    NAME Boost
    VERSION ${ogs.minimum_version.boost}
    URL https://gitlab.opengeosys.org/ogs/libs/boost-subset/-/jobs/303158/artifacts/raw/ogs-boost-${ogs.minimum_version.boost}.tar.gz
)
if(Boost_ADDED)
    add_library(Boost::boost INTERFACE IMPORTED)
    target_include_directories(Boost::boost INTERFACE "${Boost_SOURCE_DIR}")
else()
    target_include_directories(Boost::boost INTERFACE "${Boost_INCLUDE_DIR}")
endif()

CPMFindPackage(
    NAME LibXml2
    GITHUB_REPOSITORY GNOME/libxml2
    VERSION ${ogs.minimum_version.libxml2}
    OPTIONS "BUILD_SHARED_LIBS OFF"
            "LIBXML2_WITH_TESTS OFF"
            "LIBXML2_WITH_PROGRAMS OFF"
            "LIBXML2_WITH_ICONV OFF"
            "LIBXML2_WITH_ICU OFF"
            "LIBXML2_WITH_LZMA OFF"
            "LIBXML2_WITH_PYTHON OFF"
            "LIBXML2_WITH_ZLIB OFF"
    EXCLUDE_FROM_ALL YES
)
if(LibXml2_ADDED)
    add_library(LibXml2::LibXml2 ALIAS LibXml2)
    set(LIBXML2_INCLUDE_DIR ${LibXml2_SOURCE_DIR})
    list(APPEND DISABLE_WARNINGS_TARGETS LibXml2)
endif()

if(GUIX_BUILD)
    find_library(XMLPATCH_LIB xmlpatch REQUIRED)
    add_library(xmlpatch INTERFACE)
    find_package(LibXml2 REQUIRED)
    target_link_libraries(xmlpatch INTERFACE ${XMLPATCH_LIB} LibXml2::LibXml2)
else()
    CPMAddPackage(
        NAME xmlpatch
        VERSION 0.4.2
        GIT_REPOSITORY https://gitlab.opengeosys.org/ogs/libs/xmlpatch.git
        OPTIONS "BUILD_SHARED_LIBS OFF"
    )
endif()

if(OGS_BUILD_SWMM)
    CPMAddPackage(
        NAME SWMMInterface
        GITHUB_REPOSITORY ufz/SwmmInterface
        GIT_TAG 141e05ae1f419918799d7bf9178ebcd97feb1ed3
        OPTIONS "BUILD_SHARED_LIBS OFF"
    )
    if(SWMMInterface_ADDED)
        target_include_directories(
            SwmmInterface SYSTEM PUBLIC ${SWMMInterface_SOURCE_DIR}
        )
        list(APPEND DISABLE_WARNINGS_TARGETS SWMM SwmmInterface)
    endif()
endif()

CPMFindPackage(
    NAME nlohmann_json
    VERSION ${ogs.minimum_version.json}
    # the git repo is incredibly large, so we download the archived include
    # directory
    URL https://github.com/nlohmann/json/releases/download/v${ogs.minimum_version.json}/include.zip
    URL_HASH SHA256=${ogs.minimum_version.json_sha}
)
if(nlohmann_json_ADDED)
    add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED)
    target_include_directories(
        nlohmann_json::nlohmann_json
        INTERFACE ${nlohmann_json_SOURCE_DIR}/include
    )
endif()

if(OGS_BUILD_GUI)
    CPMAddPackage(
        NAME rapidxml
        VERSION 1.13
        GITHUB_REPOSITORY ufz/rapidxml
        GIT_TAG 2ae4b2888165a393dfb6382168825fddf00c27b9
        EXCLUDE_FROM_ALL YES
    )
    if(rapidxml_ADDED)
        add_library(rapidxml INTERFACE IMPORTED)
        target_include_directories(rapidxml INTERFACE ${rapidxml_SOURCE_DIR})
    endif()
endif()

if(OGS_BUILD_GUI)
    find_package(Shapelib QUIET)
    if(Shapelib_FOUND)
        add_library(shp INTERFACE IMPORTED)
        target_include_directories(shp INTERFACE ${Shapelib_INCLUDE_DIRS})
        target_link_libraries(shp INTERFACE ${Shapelib_LIBRARIES})
    else()
        CPMAddPackage(
            NAME Shapelib
            GITHUB_REPOSITORY OSGeo/shapelib
            VERSION 1.5.0-dev
            GIT_TAG 21ae8fc16afa15a1b723077b6cec3a9abc592f6a
            EXCLUDE_FROM_ALL YES
        )
        target_include_directories(
            shp INTERFACE $<BUILD_INTERFACE:${Shapelib_SOURCE_DIR}>
        )
    endif()
endif()

if(OGS_USE_CVODE)
    CPMAddPackage(
        NAME CVODE
        GITHUB_REPOSITORY ufz/cvode
        VERSION 2.8.2
        GIT_TAG 42d786bff4f950045d2de941677ecd4432cec855
        OPTIONS "EXAMPLES_ENABLE OFF"
        EXCLUDE_FROM_ALL YES
    )
    if(CVODE_ADDED)
        add_library(CVODE::CVODE INTERFACE IMPORTED)
        target_include_directories(
            CVODE::CVODE INTERFACE ${CVODE_SOURCE_DIR}/include
                                   ${CVODE_BINARY_DIR}/include
        )
        target_link_libraries(
            CVODE::CVODE INTERFACE sundials_cvode_static
                                   sundials_nvecserial_static
        )
    endif()
endif()

# VTK ###
# ~~~
# if(OGS_USE_INSITU)
#     find_package(ParaView REQUIRED)
# endif()
# ~~~
if(GUIX_BUILD)
    add_library(exprtk INTERFACE IMPORTED)
else()
    CPMAddPackage(
        NAME exprtk
        GIT_REPOSITORY https://gitlab.opengeosys.org/ogs/libs/exprtk.git
        GIT_TAG 2a5c62b93c9661470e69be572f22d821308b6f61
        DOWNLOAD_ONLY YES
    )
    if(exprtk_ADDED)
        add_library(exprtk INTERFACE IMPORTED)
        target_include_directories(exprtk SYSTEM INTERFACE ${exprtk_SOURCE_DIR})
    endif()
endif()

if(GUIX_BUILD)
    find_package(range-v3 REQUIRED)
else()
    CPMFindPackage(
        NAME range-v3
        GITHUB_REPOSITORY ericniebler/range-v3
        VERSION ${ogs.minimum_version.range-v3}
        GIT_TAG ${ogs.minimum_version.range-v3}
        EXCLUDE_FROM_ALL YES
    )
endif()

if((OGS_BUILD_TESTING OR OGS_BUILD_UTILS) AND NOT GUIX_BUILD)
    CPMAddPackage(
        NAME vtkdiff GITHUB_REPOSITORY ufz/vtkdiff
        GIT_TAG 9754b4da43c6adfb65d201ed920b5f6ea27b38b9
    )
    if(vtkdiff_ADDED)
        install(PROGRAMS $<TARGET_FILE:vtkdiff> DESTINATION bin)
    endif()
endif()

if(OGS_USE_PETSC)
    include(CheckCXXSymbolExists)
    set(CMAKE_REQUIRED_INCLUDES "${HDF5_INCLUDE_DIRS}")
    set(CMAKE_REQUIRED_LIBRARIES "${HDF5_LIBRARIES}")
    check_cxx_symbol_exists(H5Pset_fapl_mpio hdf5.h HAVE_H5Pset_fapl_mpio)
    unset(CMAKE_REQUIRED_INCLUDES)
    if(NOT HAVE_H5Pset_fapl_mpio)
        message(FATAL_ERROR "HDF5 was not build with MPI support! "
                            "(Enable with HDF5_ENABLE_PARALLEL)"
        )
    endif()
endif()

if((OGS_BUILD_TESTING OR OGS_BUILD_UTILS) AND NOT GUIX_BUILD)
    set(XDMF_LIBNAME OgsXdmf CACHE STRING "")
    CPMAddPackage(
        NAME xdmf
        VERSION 3.0.0
        GIT_REPOSITORY https://gitlab.opengeosys.org/ogs/xdmflib.git
        GIT_TAG 92a851f1acb87ad5367eb62f9b97785bedb700bb
        OPTIONS "XDMF_LIBNAME OgsXdmf" "CMAKE_MACOSX_RPATH ON"
                "HDF5_C_INCLUDE_DIR ${HDF5_INCLUDE_DIRS}"
        EXCLUDE_FROM_ALL YES
    )
    if(xdmf_ADDED)
        target_include_directories(
            OgsXdmf PUBLIC ${xdmf_SOURCE_DIR} ${xdmf_BINARY_DIR}
        )

        target_link_libraries(OgsXdmf Boost::boost)
        target_include_directories(
            OgsXdmfCore PUBLIC ${xdmf_SOURCE_DIR}/core ${xdmf_BINARY_DIR}/core
            PRIVATE ${xdmf_SOURCE_DIR}/CMake/VersionSuite
        )
        target_link_libraries(
            OgsXdmfCore PUBLIC Boost::boost LibXml2::LibXml2 ${HDF5_LIBRARIES}
        )

        set_target_properties(
            OgsXdmf OgsXdmfCore
            PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                       ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
                       LIBRARY_OUTPUT_DIRECTORY
                       ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
                       ARCHIVE_OUTPUT_DIRECTORY
                       ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
        )
        if(BUILD_SHARED_LIBS)
            install(TARGETS OgsXdmf OgsXdmfCore
                    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            )
        endif()
        list(APPEND DISABLE_WARNINGS_TARGETS OgsXdmf OgsXdmfCore)
    endif()
endif()

if(MSVC)
    CPMAddPackage(
        NAME GroupSourcesByFolder.cmake
        GITHUB_REPOSITORY TheLartians/GroupSourcesByFolder.cmake VERSION 1.0
    )
endif()

if(OGS_BUILD_UTILS)
    if(NOT GUIX_BUILD)
        set(_metis_options "MSVC ${WIN32}")
        if(WIN32)
            list(APPEND _metis_options "BUILD_SHARED_LIBS OFF")
        else()
            list(APPEND _metis_options
                 "CMAKE_C_FLAGS -D_POSIX_C_SOURCE=200809L"
            )
        endif()
        CPMFindPackage(
            NAME GKlib
            GIT_REPOSITORY https://github.com/KarypisLab/GKlib
            GIT_TAG 8bd6bad750b2b0d90800c632cf18e8ee93ad72d7
            VERSION 5.1.1
            OPTIONS "CMAKE_POLICY_DEFAULT_CMP0042 NEW" ${_metis_options}
            EXCLUDE_FROM_ALL YES
        )
        CPMFindPackage(
            NAME metis
            GIT_REPOSITORY https://github.com/KarypisLab/METIS
            VERSION 5.2.1
            EXCLUDE_FROM_ALL YES UPDATE_DISCONNECTED ON
            PATCH_COMMAND git apply
                          ${PROJECT_SOURCE_DIR}/scripts/cmake/metis.patch
            OPTIONS ${_metis_options}
        )
        if(GKlib_ADDED AND metis_ADDED)
            target_include_directories(
                metis SYSTEM
                PUBLIC ${GKlib_SOURCE_DIR} ${metis_SOURCE_DIR}/include
                       ${metis_SOURCE_DIR}/libmetis
            )
            target_compile_definitions(
                metis PUBLIC IDXTYPEWIDTH=64 REALTYPEWIDTH=32
            )
            list(APPEND DISABLE_WARNINGS_TARGETS metis GKlib)
            install(TARGETS mpmetis GKlib)
        endif()
    else()
        find_library(METIS_LIB metis REQUIRED)
        find_path(METIS_INC "metis.h" REQUIRED)
        find_program(MPMETIS_TOOL mpmetis REQUIRED)
        message(STATUS "Metis: ${METIS_LIB} | ${METIS_INC} | ${MPMETIS_TOOL}")
        add_library(metis INTERFACE IMPORTED)
        target_include_directories(metis SYSTEM INTERFACE ${METIS_INC})
        target_link_libraries(metis INTERFACE ${METIS_LIB})
    endif()
endif()

if(OGS_USE_NETCDF)
    find_package(netCDF CONFIG REQUIRED)
    find_library(NETCDF_LIBRARIES_CXX NAMES netcdf_c++4 netcdf-cxx4)
    if(NOT NETCDF_LIBRARIES_CXX)
        CPMAddPackage(
            NAME netcdf-cxx4
            GIT_REPOSITORY https://github.com/Unidata/netcdf-cxx4
            VERSION 4.3.1
            EXCLUDE_FROM_ALL YES SOURCE_SUBDIR cxx4
            OPTIONS "NCXX_ENABLE_TESTS OFF"
        )
        set_target_properties(
            netCDF::netcdf PROPERTIES INTERFACE_LINK_LIBRARIES ""
        ) # fix win installed config
        target_link_libraries(netcdf-cxx4 netCDF::netcdf)
    else()
        find_path(NETCDF_INCLUDES_CXX NAMES netcdf)
        add_library(netcdf-cxx4 INTERFACE IMPORTED)
        target_include_directories(netcdf-cxx4 INTERFACE ${NETCDF_INCLUDES_CXX})
        target_link_libraries(
            netcdf-cxx4 INTERFACE ${NETCDF_LIBRARIES_CXX} netCDF::netcdf
        )
    endif()
endif()

# Disable warnings
if(WIN32 AND VTK_ADDED)
    list(APPEND DISABLE_WARNINGS_TARGETS vtksys)
endif()
foreach(target ${DISABLE_WARNINGS_TARGETS})
    target_compile_options(
        ${target} PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
                          $<$<CXX_COMPILER_ID:MSVC>:/W0>
    )
endforeach()

# Hack: Disable tests from dependencies
configure_file(
    ${PROJECT_SOURCE_DIR}/scripts/cmake/test/CTestCustom.in.cmake
    ${PROJECT_BINARY_DIR}/CTestCustom.cmake @ONLY
)

find_program(CLANG_FORMAT_PROGRAM clang-format)
find_program(CMAKE_FORMAT_PROGRAM cmake-format)

if(CLANG_FORMAT_PROGRAM OR CMAKE_FORMAT_PROGRAM)
    if(NOT CMAKE_FORMAT_PROGRAM)
        set(_skip_cmake "FORMAT_SKIP_CMAKE YES")
    endif()
    CPMAddPackage(
        NAME Format.cmake
        VERSION 1.7.0
        GITHUB_REPOSITORY TheLartians/Format.cmake
        OPTIONS
            ${_skip_cmake}
            "CMAKE_FORMAT_EXCLUDE scripts/cmake/CPM.cmake|.*/Tests.cmake|scripts/cmake/vector-of-bool/.*"
    )
endif()

# Third-party licenses
set(_licenses_file ${PROJECT_BINARY_DIR}/third_party_licenses.txt)
if(NOT EXISTS ${_licenses_file})
    set(_licenses_string "")
    # Adapted from https://github.com/cpm-cmake/CPMLicenses.cmake:
    set(_print_delimiter OFF)
    foreach(package ${CPM_PACKAGES} ${_EXT_LIBS})
        file(
            GLOB
            licenses
            "${${package}_SOURCE_DIR}/LICENSE*"
            "${${package}_SOURCE_DIR}/LICENCE*"
            "${${package}_SOURCE_DIR}/COPYING*"
            "${${package}_SOURCE_DIR}/Copyright*"
        )
        list(LENGTH licenses LICENSE_COUNT)
        if(LICENSE_COUNT GREATER_EQUAL 1)
            if(_print_delimiter)
                set(_licenses_string
                    "${_licenses_string}\n----------------------------------------------------------------------------\n\n"
                )
            endif()

            list(GET licenses 0 _license)
            file(READ ${_license} license_TEXT)
            set(_licenses_string
                "${_licenses_string}The following software may be included in this product: **${package}**.\nThis software contains the following license and notice below:\n\n${license_TEXT}\n"
            )
            set(_print_delimiter ON)
        else()
            message(
                VERBOSE
                "WARNING: no license files found for package \"${package}\" in ${${package}_SOURCE_DIR} ."
            )
        endif()
        if(LICENSE_COUNT GREATER 1)
            message(
                VERBOSE
                "WARNING: multiple license files found for package \"${package}\": ${licenses}. Only first file will be used."
            )
        endif()
    endforeach()

    file(READ
         ${PROJECT_SOURCE_DIR}/scripts/cmake/packaging/additional-licenses.txt
         _additional_licenses
    )

    file(WRITE ${_licenses_file} "${_licenses_string}\n${_additional_licenses}")
    message(STATUS "Wrote licenses to ${_licenses_file}")
endif()

unset(CMAKE_FOLDER)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "└─ End Dependencies.cmake")
