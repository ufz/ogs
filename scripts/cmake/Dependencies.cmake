message(STATUS "┌─ Dependencies.cmake")
list(APPEND CMAKE_MESSAGE_INDENT "│    ")
set(CMAKE_FOLDER ThirdParty)

# ccache, on Windows requires https://github.com/cristianadam/ccache/releases
if(NOT WIN32 AND CCACHE_TOOL_PATH AND NOT OGS_DISABLE_COMPILER_CACHE)
    set(CCACHE_OPTIONS "CCACHE_SLOPPINESS=pch_defines,time_macros")
    if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang|AppleClang")
        list(APPEND CCACHE_OPTIONS "CCACHE_CPP2=true")
    endif()
    CPMAddPackage(
        NAME Ccache.cmake
        GITHUB_REPOSITORY TheLartians/Ccache.cmake
        VERSION 1.2.2
        OPTIONS "USE_CCACHE ON"
    )
elseif(
    WIN32
    AND CCACHE_TOOL_PATH
    AND NOT OGS_DISABLE_COMPILER_CACHE
    AND "${CMAKE_GENERATOR}" STREQUAL "Ninja"
)
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_TOOL_PATH} CACHE STRING "" FORCE)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_TOOL_PATH} CACHE STRING "" FORCE)
    message(STATUS "Using ccache (${CCACHE_TOOL_PATH}).")
endif()

if(OGS_BUILD_TESTING)
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        VERSION ${ogs.minimum_version.gtest}
        GIT_TAG release-${ogs.minimum_version.gtest}
        OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt ON"
                "BUILD_SHARED_LIBS OFF"
        EXCLUDE_FROM_ALL YES
    )

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

CPMFindPackage(NAME spdlog GITHUB_REPOSITORY gabime/spdlog VERSION 1.8.2)

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

CPMAddPackage(
    NAME tetgen GITHUB_REPOSITORY ufz/tetgen
    GIT_TAG 603ba181ebfaed38eec88532720e282606009b73
)
if(tetgen_ADDED)
    install(PROGRAMS $<TARGET_FILE:tetgen> DESTINATION bin)
    list(APPEND DISABLE_WARNINGS_TARGETS tet tetgen)
endif()

if(OGS_USE_PYTHON)
    CPMAddPackage(
        NAME pybind11 GITHUB_REPOSITORY pybind/pybind11
        GIT_TAG f1abf5d9159b805674197f6bc443592e631c9130
    )
endif()

if(_build_chemistry_lib)
    CPMAddPackage(
        NAME iphreeqc
        GITHUB_REPOSITORY ufz/iphreeqc
        GIT_TAG b1047d3eb03e7ef1b850231be35acb9c6a2cf345
        DOWNLOAD_ONLY YES
    )
    if(iphreeqc_ADDED)
        include(scripts/cmake/iphreeqc.cmake)
        list(APPEND DISABLE_WARNINGS_TARGETS iphreeqc)
    endif()
endif()

CPMFindPackage(
    NAME Eigen3
    # Error as in
    # https://gitlab.com/gitlab-com/gl-infra/reliability/-/issues/8475
    # GITLAB_REPOSITORY libeigen/eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/${ogs.minimum_version.eigen}/eigen-${ogs.minimum_version.eigen}.tar.gz
    GIT_TAG ${ogs.minimum_version.eigen}
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
    CPMAddPackage(
        NAME MGIS
        GITHUB_REPOSITORY ufz/MFrontGenericInterfaceSupport
        GIT_TAG ff51a694957b2f1c16689962529fad0ee2b04558
        OPTIONS "enable-doxygen-doc OFF" "enable-fortran-bindings OFF"
        EXCLUDE_FROM_ALL YES
    )
    if(MGIS_ADDED)
        set_target_properties(MFrontGenericInterface PROPERTIES CXX_STANDARD 11)
        list(APPEND DISABLE_WARNINGS_TARGETS MFrontGenericInterface)
        set(_MFRONT_TFEL_FOUND ON CACHE INTERNAL "")
    endif()
endif()

string(REPLACE "." "_" BOOST_VERSION_UNDESCORE ${ogs.minimum_version.boost})
CPMFindPackage(
    NAME Boost
    VERSION ${ogs.minimum_version.boost}
    URL https://boostorg.jfrog.io/artifactory/main/release/${ogs.minimum_version.boost}/source/boost_${BOOST_VERSION_UNDESCORE}.tar.gz
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

CPMAddPackage(
    NAME xmlpatch
    VERSION 0.4.2
    GIT_REPOSITORY https://gitlab.opengeosys.org/ogs/libs/xmlpatch.git
    OPTIONS "BUILD_SHARED_LIBS OFF"
)

if(OGS_BUILD_SWMM)
    CPMAddPackage(
        NAME SWMMInterface GITHUB_REPOSITORY ufz/SwmmInterface
        GIT_TAG 141e05ae1f419918799d7bf9178ebcd97feb1ed3
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
    VERSION 3.6.1
    # the git repo is incredibly large, so we download the archived include
    # directory
    URL https://github.com/nlohmann/json/releases/download/v3.6.1/include.zip
    URL_HASH
        SHA256=69cc88207ce91347ea530b227ff0776db82dcb8de6704e1a3d74f4841bc651cf
)
if(nlohmann_json_ADDED)
    add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED)
    target_include_directories(
        nlohmann_json::nlohmann_json INTERFACE ${nlohmann_json_SOURCE_DIR}
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
unset(VTK_OPTIONS)
foreach(option_index ${ogs.libraries.vtk.options})
    if(${ogs.libraries.vtk.options_${option_index}.condition.cmake})
        foreach(cmake_index ${ogs.libraries.vtk.options_${option_index}.cmake})
            string(
                REPLACE
                    "="
                    " "
                    cmake_option
                    "${ogs.libraries.vtk.options_${option_index}.cmake_${cmake_index}}"
            )
            list(APPEND VTK_OPTIONS ${cmake_option})
        endforeach()
    endif()
endforeach()
list(REMOVE_DUPLICATES VTK_OPTIONS)

# TODO: if(OGS_INSITU) find_package(ParaView REQUIRED) end()
unset(VTK_COMPONENTS)
foreach(opt ${VTK_OPTIONS})
    if("${opt}" MATCHES "^VTK_MODULE_ENABLE_VTK_(.*) YES")
        list(APPEND VTK_COMPONENTS ${CMAKE_MATCH_1})
    endif()
endforeach()
find_package(VTK ${ogs.minimum_version.vtk} QUIET COMPONENTS ${VTK_COMPONENTS})

if(NOT VTK_FOUND)
    list(APPEND VTK_OPTIONS "BUILD_SHARED_LIBS OFF")
    # Workaround for configuration error in [vtk]/CMake/vtkGroups.cmake:43
    set(VTK_GROUP_ENABLE_Rendering DONT_WANT CACHE STRING "")
    set(VTK_GROUP_ENABLE_StandAlone DONT_WANT CACHE STRING "")

    CPMAddPackage(
        NAME VTK
        GITHUB_REPOSITORY kitware/vtk
        VERSION ${ogs.minimum_version.vtk}
        OPTIONS ${VTK_OPTIONS}
        EXCLUDE_FROM_ALL YES GIT_SUBMODULES "" # Disable submodules
    )
endif()
if(VTK_ADDED)
    if(OpenMP_FOUND AND TARGET vtkFiltersStatistics)
        target_link_libraries(vtkFiltersStatistics PRIVATE OpenMP::OpenMP_C)
    endif()
    if(TARGET loguru)
        # Fixes https://stackoverflow.com/questions/9894961 on vismac05:
        set_target_properties(loguru PROPERTIES CXX_VISIBILITY_PRESET default)
        # Also suppress warnings
        list(APPEND DISABLE_WARNINGS_TARGETS loguru)
    endif()
endif()
# end VTK ###

if(VTK_ADDED)
    # VTK already comes with exprtk, reusing it.
    target_include_directories(
        exprtk SYSTEM
        INTERFACE
            $<BUILD_INTERFACE:${VTK_SOURCE_DIR}/ThirdParty/exprtk/vtkexprtk>
    )
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

if(OGS_BUILD_TESTING OR OGS_BUILD_UTILS)
    CPMAddPackage(
        NAME vtkdiff GITHUB_REPOSITORY ufz/vtkdiff
        GIT_TAG 788100291f73e472febf7e5550eea36ec4be518b
    )
    if(vtkdiff_ADDED)
        install(PROGRAMS $<TARGET_FILE:vtkdiff> DESTINATION bin)
    endif()
endif()

if(OGS_USE_MPI)
    set(_hdf5_options "HDF5_ENABLE_PARALLEL ON")
endif()

string(REPLACE "." "_" HDF5_TAG ${ogs.minimum_version.hdf5})
if(OGS_USE_NETCDF)
    list(APPEND CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR})
    find_package(HDF5 REQUIRED)
else()
    # ZLIB is a HDF5 dependency
    if(NOT VTK_ADDED)
        CPMFindPackage(
            NAME ZLIB
            GITHUB_REPOSITORY madler/zlib
            VERSION 1.2.11
            EXCLUDE_FROM_ALL YES
        )
    endif()

    CPMFindPackage(
        NAME HDF5
        GITHUB_REPOSITORY HDFGroup/hdf5
        GIT_TAG hdf5-${HDF5_TAG}
        VERSION ${ogs.minimum_version.hdf5}
        OPTIONS "HDF5_EXTERNALLY_CONFIGURED 1"
                "HDF5_GENERATE_HEADERS OFF"
                "HDF5_BUILD_TOOLS OFF"
                "HDF5_BUILD_EXAMPLES OFF"
                "HDF5_BUILD_HL_LIB OFF"
                "HDF5_BUILD_FORTRAN OFF"
                "HDF5_BUILD_CPP_LIB OFF"
                "HDF5_BUILD_JAVA OFF"
                ${_hdf5_options}
        EXCLUDE_FROM_ALL YES
    )
    if(HDF5_ADDED)
        list(APPEND DISABLE_WARNINGS_TARGETS hdf5-static)
        set(HDF5_LIBRARIES hdf5-static)
        if(ZLIB_ADDED)
            list(APPEND HDF5_LIBRARIES zlibstatic)
        endif()
        set(HDF5_INCLUDE_DIRS ${HDF5_SOURCE_DIR}/src ${HDF5_BINARY_DIR})
        set(HDF5_C_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
        set(HDF5_C_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})
        target_include_directories(hdf5-static INTERFACE ${HDF5_INCLUDE_DIRS})
        if(OGS_USE_MKL AND WIN32)
            # In this case the hdf5 build fails with, e.g.:
            # ~~~
            # H5system.c(710): error C2065: 'timezone': undeclared identifier
            # ~~~
            # Reason is that the H5_HAVE_VISUAL_STUDIO-symbol is not defined
            # anymore!?
            target_compile_definitions(
                hdf5-static PRIVATE -DH5_HAVE_VISUAL_STUDIO
            )
        endif()
    else()
        find_package(HDF5 REQUIRED)
    endif()
endif()

if(OGS_USE_PETSC AND NOT HDF5_ADDED)
    include(CheckCXXSymbolExists)
    set(CMAKE_REQUIRED_INCLUDES "${HDF5_INCLUDE_DIR}" "${HDF5_BINARY_DIR}")
    set(CMAKE_REQUIRED_LIBRARIES "${HDF5_LIBRARIES}")
    check_cxx_symbol_exists(H5Pset_fapl_mpio hdf5.h HAVE_H5Pset_fapl_mpio)
    unset(CMAKE_REQUIRED_INCLUDES)
    if(NOT HAVE_H5Pset_fapl_mpio)
        message(FATAL_ERROR "HDF5 was not build with MPI support! "
                            "(Enable with HDF5_ENABLE_PARALLEL)"
        )
    endif()
endif()

# Does not compile in Debug-mode, see #3175.
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND OGS_BUILD_TESTING)
    set(XDMF_LIBNAME OgsXdmf CACHE STRING "")
    CPMAddPackage(
        NAME xdmf
        VERSION 3.0.0
        GIT_REPOSITORY https://gitlab.opengeosys.org/ogs/xdmflib.git
        GIT_TAG 92a851f1acb87ad5367eb62f9b97785bedb700bb
        OPTIONS "XDMF_LIBNAME OgsXdmf"
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

if(OGS_BUILD_UTILS)
    CPMAddPackage(
        NAME metis
        GITHUB_REPOSITORY scibuilder/metis
        GIT_TAG 982842a5ace9b3da2b2800817eb9e5fd3b42966b
        DOWNLOAD_ONLY YES
    )
    include(${PROJECT_SOURCE_DIR}/scripts/cmake/MetisSetup.cmake)
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
            "CMAKE_FORMAT_EXCLUDE scripts/cmake/CPM.cmake|.*/Tests.cmake|scripts/cmake/jedbrown/.*|scripts/cmake/conan/conan.cmake|scripts/cmake/vector-of-bool/.*"
    )
endif()

# Third-party licenses
CPMAddPackage(
    NAME CPMLicenses.cmake GITHUB_REPOSITORY cpm-cmake/CPMLicenses.cmake
    VERSION 0.0.5
)
cpm_licenses_create_disclaimer_target(
    write-licenses "${PROJECT_BINARY_DIR}/third_party_licenses.txt"
    "${CPM_PACKAGES}"
)

CPMAddPackage(
    NAME GroupSourcesByFolder.cmake
    GITHUB_REPOSITORY TheLartians/GroupSourcesByFolder.cmake VERSION 1.0
)

unset(CMAKE_FOLDER)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "└─ End Dependencies.cmake")
