message(STATUS "┌─ Dependencies.cmake")
list(APPEND CMAKE_MESSAGE_INDENT "│    ")
set(CMAKE_FOLDER ThirdParty)

if(OGS_BUILD_TESTING)
    if(GUIX_BUILD)
        find_package(GTest REQUIRED)
        add_library(gtest ALIAS GTest::gtest)
        add_library(gmock ALIAS GTest::gmock)
        find_program(VTKDIFF_TOOL vtkdiff REQUIRED)
        add_executable(vtkdiff IMPORTED GLOBAL)
        set_target_properties(vtkdiff PROPERTIES
            IMPORTED_LOCATION "${VTKDIFF_TOOL}"
        )
        add_library(autocheck INTERFACE IMPORTED)
    else()
        CPMAddPackage(
            NAME googletest
            GITHUB_REPOSITORY google/googletest
            VERSION ${ogs.minimum_version.gtest}
            GIT_TAG v${ogs.tested_version.gtest}
            OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt ON"
                    "BUILD_SHARED_LIBS OFF"
            EXCLUDE_FROM_ALL YES SYSTEM TRUE
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
endif()

if(GUIX_BUILD)
    find_package(spdlog REQUIRED)
    target_compile_definitions(spdlog::spdlog INTERFACE SPDLOG_FMT_EXTERNAL)
else()
    CPMAddPackage(
        NAME spdlog
        GITHUB_REPOSITORY gabime/spdlog
        VERSION 1.15.3
        OPTIONS "BUILD_SHARED_LIBS OFF" SYSTEM TRUE
    )
endif()

if(GUIX_BUILD OR CONDA_BUILD)
    find_path(_tclap_include UnlabeledValueArg.h PATH_SUFFIXES tclap REQUIRED)
    add_library(tclap INTERFACE IMPORTED)
    target_include_directories(
        tclap SYSTEM INTERFACE ${_tclap_include} ${_tclap_include}/..
    )
else()
    CPMFindPackage(
        NAME tclap
        GIT_REPOSITORY https://git.code.sf.net/p/tclap/code
        GIT_TAG 81b3d2a0c47895c22e9bb8c577f5ab521f76e5d2
        VERSION 1.4.0
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
    find_package(pybind11 REQUIRED)
else()
    set(PYBIND11_FINDPYTHON ON)
    CPMFindPackage(
        NAME pybind11 GITHUB_REPOSITORY pybind/pybind11
        VERSION ${ogs.minimum_version.pybind11} SYSTEM TRUE
        OPTIONS "CMAKE_POLICY_VERSION_MINIMUM 3.10"
    )
endif()

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
            OPTIONS "CMAKE_POLICY_VERSION_MINIMUM 3.10"
        )
        if(iphreeqc_ADDED)
            target_include_directories(
                IPhreeqc SYSTEM
                PUBLIC ${iphreeqc_SOURCE_DIR}/src
                INTERFACE ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/common
                          ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/PhreeqcKeywords
                          ${iphreeqc_SOURCE_DIR}/src/phreeqcpp
            )
            if(BUILD_SHARED_LIBS)
                install(TARGETS IPhreeqc
                        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                )
            endif()
        endif()
        list(APPEND DISABLE_WARNINGS_TARGETS IPhreeqc)
    endif()
endif()

if(GUIX_BUILD)
    find_package(Eigen3 REQUIRED)
    if("${Eigen3_VERSION}" VERSION_GREATER 3.4)
        set(OGS_USE_EIGEN_UNSUPPORTED ON CACHE BOOL "" FORCE)
    else()
        set(OGS_USE_EIGEN_UNSUPPORTED OFF CACHE BOOL "" FORCE)
    endif()
endif()

set(_eigen_version ${ogs.minimum_version.eigen})
set(_eigen_url
    https://gitlab.com/libeigen/eigen/-/archive/${_eigen_version}/eigen-${_eigen_version}.tar.gz
)
if(OGS_EIGEN_PARALLEL_BACKEND STREQUAL "MKL" AND NOT OGS_USE_EIGEN_UNSUPPORTED)
    message(
        FATAL_ERROR
            "OGS_EIGEN_PARALLEL_BACKEND=MKL requires OGS_USE_EIGEN_UNSUPPORTED!"
    )
endif()

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
    )
    if(${OGS_EIGEN_PARALLEL_BACKEND} STREQUAL "OpenMP")
        target_include_directories(
            Eigen3::Eigen SYSTEM INTERFACE ${OpenMP_CXX_INCLUDE_DIRS}
        )
    else()
        target_include_directories(
            Eigen3::Eigen SYSTEM INTERFACE ${MKL_INCLUDE}
        )
    endif()
endif()

if(OGS_USE_MFRONT)
    if(GUIX_BUILD OR CONDA_BUILD)
        if(CONDA_BUILD)
            set(MGIS_DIR $ENV{CONDA_PREFIX})
        endif()
        find_package(MFrontGenericInterface REQUIRED)
        add_library(MFrontGenericInterface ALIAS mgis::MFrontGenericInterface)
    else()
        set(CMAKE_REQUIRE_FIND_PACKAGE_TFEL TRUE)
        # Patch only works when CPM_SOURCE_CACHE is set. Following conditional
        # logic can be removed if
        # https://github.com/cpm-cmake/CPM.cmake/issues/577 is resolved.
        if(NOT "${CPM_SOURCE_CACHE}" STREQUAL "")
            set(_mgis_patch_args
                PATCHES ${PROJECT_SOURCE_DIR}/scripts/cmake/mgis-flags.patch
            )
            message(STATUS "Adding mgis-flags.patch.")
        endif()
        if(OGS_CPU_ARCHITECTURE STREQUAL "generic")
            set(_enable_portable_build_option "enable-portable-build ON")
        endif()
        CPMAddPackage(
            NAME MGIS
            GITHUB_REPOSITORY thelfer/MFrontGenericInterfaceSupport
            GIT_TAG rliv-2.2
            OPTIONS "enable-doxygen-doc OFF" "enable-fortran-bindings OFF"
                    "enable-website OFF" "CMAKE_POLICY_VERSION_MINIMUM 3.10" ${_enable_portable_build_option}
            EXCLUDE_FROM_ALL YES SYSTEM TRUE ${_mgis_patch_args}
        )
        list(APPEND DISABLE_WARNINGS_TARGETS MFrontGenericInterface)
    endif()
endif()

# Boost libraries used by ogs, can be linked with Boost::[lib_name]
set(BOOST_INCLUDE_LIBRARIES
    math
    property_tree
    algorithm
    smart_ptr
    tokenizer
    assign
    dynamic_bitset
    range
    variant
    interprocess
)
if(GUIX_BUILD)
    find_package(Boost REQUIRED)
else()
    if(OGS_BUILD_WHEEL)
        set(_boost_options "BUILD_SHARED_LIBS OFF")
    endif()
    CPMFindPackage(
        NAME Boost
        VERSION ${ogs.minimum_version.boost}
        URL https://github.com/boostorg/boost/releases/download/boost-${ogs.minimum_version.boost}/boost-${ogs.minimum_version.boost}.tar.xz
        OPTIONS "BOOST_ENABLE_CMAKE ON"
                "CMAKE_POLICY_VERSION_MINIMUM 3.10"
                ${_boost_options}
    )
endif()
if(NOT Boost_ADDED)
    # Boost from system found. There are only Boost::headers and Boost::boost
    # targets.
    foreach(lib ${BOOST_INCLUDE_LIBRARIES})
        add_library(Boost::${lib} ALIAS Boost::headers)
    endforeach()
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
    EXCLUDE_FROM_ALL YES SYSTEM TRUE
)
if(LibXml2_ADDED)
    add_library(LibXml2::LibXml2 ALIAS LibXml2)
    set(LIBXML2_INCLUDE_DIR ${LibXml2_SOURCE_DIR})
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
        OPTIONS "BUILD_SHARED_LIBS OFF" SYSTEM TRUE
    )
endif()

if(OGS_BUILD_SWMM)
    CPMAddPackage(
        NAME SWMMInterface
        GITHUB_REPOSITORY ufz/SwmmInterface
        GIT_TAG 141e05ae1f419918799d7bf9178ebcd97feb1ed3
        OPTIONS "BUILD_SHARED_LIBS OFF" SYSTEM TRUE
    )
    if(SWMMInterface_ADDED)
        target_include_directories(
            SwmmInterface SYSTEM PUBLIC ${SWMMInterface_SOURCE_DIR}
        )
    endif()
    list(APPEND DISABLE_WARNINGS_TARGETS SWMM SwmmInterface)
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
        nlohmann_json::nlohmann_json SYSTEM
        INTERFACE ${nlohmann_json_SOURCE_DIR}/include
    )
endif()

if(OGS_BUILD_GUI)
    CPMAddPackage(
        NAME rapidxml
        VERSION 1.13
        GITHUB_REPOSITORY ufz/rapidxml
        GIT_TAG 2ae4b2888165a393dfb6382168825fddf00c27b9
        EXCLUDE_FROM_ALL YES SYSTEM TRUE
    )
    if(rapidxml_ADDED)
        add_library(rapidxml INTERFACE IMPORTED)
        target_include_directories(
            rapidxml SYSTEM INTERFACE ${rapidxml_SOURCE_DIR}
        )
    endif()
endif()

if(OGS_BUILD_GUI)
    find_package(Shapelib QUIET)
    if(Shapelib_FOUND)
        add_library(shp INTERFACE IMPORTED)
        target_include_directories(
            shp SYSTEM INTERFACE ${Shapelib_INCLUDE_DIRS}
        )
        target_link_libraries(shp INTERFACE ${Shapelib_LIBRARIES})
    else()
        CPMAddPackage(
            NAME Shapelib
            GITHUB_REPOSITORY OSGeo/shapelib
            VERSION 1.5.0-dev
            GIT_TAG 21ae8fc16afa15a1b723077b6cec3a9abc592f6a
            EXCLUDE_FROM_ALL YES SYSTEM TRUE
        )
        target_include_directories(
            shp SYSTEM INTERFACE $<BUILD_INTERFACE:${Shapelib_SOURCE_DIR}>
        )
    endif()
endif()

# VTK ###
# ~~~
# if(OGS_USE_INSITU)
#     find_package(ParaView REQUIRED)
# endif()
# ~~~
if(GUIX_BUILD OR CONDA_BUILD)
    add_library(exprtk INTERFACE IMPORTED)
else()
    CPMAddPackage(
        NAME exprtk
        GITHUB_REPOSITORY ArashPartow/exprtk
        VERSION 0.0.3
        GIT_TAG 0.0.3
        DOWNLOAD_ONLY YES
    )
    if(exprtk_ADDED)
        add_library(exprtk INTERFACE IMPORTED)
        target_include_directories(exprtk SYSTEM INTERFACE ${exprtk_SOURCE_DIR})
    endif()
endif()

if(GUIX_BUILD OR CONDA_BUILD)
    find_package(range-v3 REQUIRED)
else()
    CPMFindPackage(
        NAME range-v3
        GITHUB_REPOSITORY ericniebler/range-v3
        VERSION ${ogs.minimum_version.range-v3}
        GIT_TAG ${ogs.minimum_version.range-v3}
        EXCLUDE_FROM_ALL YES SYSTEM TRUE
        OPTIONS "CMAKE_POLICY_VERSION_MINIMUM 3.10"
    )
endif()

if(NOT (GUIX_BUILD OR CONDA_BUILD))
    if((OGS_BUILD_TESTING OR OGS_BUILD_UTILS))
        CPMAddPackage(
            NAME vtkdiff GITHUB_REPOSITORY ufz/vtkdiff
            GIT_TAG 628c4694783f865d7f0ab3ba9bdd5530ce4567e9
        )
        if(vtkdiff_ADDED)
            install(PROGRAMS $<TARGET_FILE:vtkdiff> DESTINATION bin)
        endif()
    endif()
endif()

if(OGS_USE_PETSC)
    include(CheckCXXSymbolExists)
    set(CMAKE_REQUIRED_INCLUDES "${HDF5_INCLUDE_DIRS}")
    set(CMAKE_REQUIRED_LIBRARIES "${HDF5_LIBRARIES}")
    check_cxx_symbol_exists(H5Pset_fapl_mpio hdf5.h HAVE_H5Pset_fapl_mpio)
    unset(CMAKE_REQUIRED_INCLUDES)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(NOT HAVE_H5Pset_fapl_mpio)
        message(FATAL_ERROR "HDF5 was not build with MPI support! "
                            "(Enable with HDF5_ENABLE_PARALLEL)"
        )
    endif()
endif()

if(NOT GUIX_BUILD)
    set(XDMF_LIBNAME OgsXdmf CACHE STRING "")
    CPMAddPackage(
        NAME xdmf
        VERSION 3.0.0
        GIT_REPOSITORY https://gitlab.opengeosys.org/ogs/xdmflib.git
        GIT_TAG 374ee63abf605ab4c6639989bebc5096881f4f57
        OPTIONS "XDMF_LIBNAME OgsXdmf" "CMAKE_MACOSX_RPATH ON"
                "HDF5_C_INCLUDE_DIR ${HDF5_INCLUDE_DIRS}"
                "CMAKE_POLICY_VERSION_MINIMUM 3.10"
        EXCLUDE_FROM_ALL YES SYSTEM TRUE
    )
    if(xdmf_ADDED)
        target_include_directories(
            OgsXdmf SYSTEM PUBLIC ${xdmf_SOURCE_DIR} ${xdmf_BINARY_DIR}
        )

        target_link_libraries(OgsXdmf Boost::tokenizer)
        target_include_directories(
            OgsXdmfCore SYSTEM PUBLIC ${xdmf_SOURCE_DIR}/core
                                      ${xdmf_BINARY_DIR}/core
            PRIVATE ${xdmf_SOURCE_DIR}/CMake/VersionSuite
        )
        target_link_libraries(
            OgsXdmfCore PUBLIC LibXml2::LibXml2 ${HDF5_LIBRARIES} Boost::variant
                               Boost::smart_ptr
            PRIVATE Boost::tokenizer Boost::assign Boost::algorithm
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
    endif()
    list(APPEND DISABLE_WARNINGS_TARGETS OgsXdmf OgsXdmfCore)
endif()

if(MSVC)
    CPMAddPackage(
        NAME GroupSourcesByFolder.cmake
        GITHUB_REPOSITORY TheLartians/GroupSourcesByFolder.cmake VERSION 1.0
    )
endif()

if(OGS_BUILD_UTILS)
    if(NOT GUIX_BUILD AND NOT CONDA_BUILD)
        set(_metis_options "MSVC ${WIN32}")
        if(WIN32)
            list(APPEND _metis_options "BUILD_SHARED_LIBS OFF")
        else()
            list(APPEND _metis_options
                 "CMAKE_C_FLAGS -D_POSIX_C_SOURCE=200809L" ${CPU_FLAGS}
            )
        endif()
        if(OGS_CPU_ARCHITECTURE STREQUAL "generic")
            set(_gklib_patch_args
                PATCHES ${PROJECT_SOURCE_DIR}/scripts/cmake/gklib.patch
            )
            set(_metis_genric_patch
                ${PROJECT_SOURCE_DIR}/scripts/cmake/metis-generic-build.patch
            )
            message(STATUS "Adding gklib.patch.")
            message(STATUS "Adding metis-generic-build.patch.")
        endif()
        CPMFindPackage(
            NAME GKlib
            GIT_REPOSITORY https://github.com/KarypisLab/GKlib
            GIT_TAG 8bd6bad750b2b0d90800c632cf18e8ee93ad72d7
            VERSION 5.1.1
            OPTIONS "CMAKE_POLICY_DEFAULT_CMP0042 NEW" ${_metis_options}
                    "CMAKE_POLICY_VERSION_MINIMUM 3.10"
            EXCLUDE_FROM_ALL YES SYSTEM TRUE ${_gklib_patch_args}
        )
        CPMFindPackage(
            NAME metis
            GIT_REPOSITORY https://github.com/KarypisLab/METIS
            VERSION 5.2.1
            EXCLUDE_FROM_ALL YES UPDATE_DISCONNECTED ON
            PATCHES ${PROJECT_SOURCE_DIR}/scripts/cmake/metis.patch ${_metis_genric_patch}
            OPTIONS ${_metis_options} "CMAKE_POLICY_VERSION_MINIMUM 3.10"
            SYSTEM TRUE
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
            install(TARGETS mpmetis GKlib)
        endif()
        list(APPEND DISABLE_WARNINGS_TARGETS metis mpmetis GKlib)
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
    if(NOT GUIX_BUILD AND NOT CONDA_BUILD)
        find_package(netCDF CONFIG REQUIRED)
        find_library(NETCDF_LIBRARIES_CXX NAMES netcdf_c++4 netcdf-cxx4)
        if(NOT NETCDF_LIBRARIES_CXX)
            CPMAddPackage(
                NAME netcdf-cxx4
                GIT_REPOSITORY https://github.com/Unidata/netcdf-cxx4
                VERSION 4.3.1
                EXCLUDE_FROM_ALL YES SOURCE_SUBDIR cxx4
                OPTIONS "NCXX_ENABLE_TESTS OFF" SYSTEM TRUE
            )
            set_target_properties(
                netCDF::netcdf PROPERTIES INTERFACE_LINK_LIBRARIES ""
            ) # fix win installed config
            target_link_libraries(netcdf-cxx4 netCDF::netcdf)
        else()
            find_path(NETCDF_INCLUDES_CXX NAMES netcdf)
            add_library(netcdf-cxx4 INTERFACE IMPORTED)
            target_include_directories(
                netcdf-cxx4 SYSTEM INTERFACE ${NETCDF_INCLUDES_CXX}
            )
            target_link_libraries(
                netcdf-cxx4 INTERFACE ${NETCDF_LIBRARIES_CXX} netCDF::netcdf
            )
        endif()
    else()
        find_library(
            NETCDF_LIBRARIES_CXX NAMES netcdf_c++4 netcdf-cxx4 REQUIRED
        )
        find_path(NETCDF_INCLUDES_CXX "ncVar.h" REQUIRED)
        message(
            STATUS
                "NetCDF-cxx4: ${NETCDF_LIBRARIES_CXX} | ${NETCDF_INCLUDES_CXX}"
        )
        add_library(netcdf-cxx4 INTERFACE IMPORTED)
        target_include_directories(
            netcdf-cxx4 SYSTEM INTERFACE ${NETCDF_INCLUDES_CXX}
        )
        target_link_libraries(netcdf-cxx4 INTERFACE ${NETCDF_LIBRARIES_CXX})
    endif()
endif()

# Disable warnings
if(WIN32 AND VTK_ADDED)
    list(APPEND DISABLE_WARNINGS_TARGETS vtksys)
endif()
foreach(target ${DISABLE_WARNINGS_TARGETS})
    if(TARGET ${target})
        target_compile_options(
            ${target} PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
                              $<$<CXX_COMPILER_ID:MSVC>:/W0>
        )
    endif()
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
