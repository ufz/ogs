# Build dependencies via ExternalProject_Add() at configure time in
# ${PROJECT_BINARY_DIR}/_ext
include(BuildExternalProject)

set(OGS_EXTERNAL_DEPENDENCIES_CACHE ""
    CACHE PATH "Directory containing source archives of external dependencies."
)

if(OGS_USE_MFRONT)
    option(OGS_BUILD_TFEL
           "Build TFEL locally. Needs to be set with a clean cache!" OFF
    )
    set(_tfel_source
        GIT_REPOSITORY
        https://github.com/${ogs.minimum_version.tfel-repo}/tfel.git GIT_TAG
        rliv-${ogs.minimum_version.tfel-rliv}
    )
    set(_tfel_source_file
        ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/tfel-rliv-${ogs.minimum_version.tfel-rliv}.zip
    )
    if(EXISTS ${_tfel_source_file})
        set(_tfel_source URL ${_tfel_source_file})
    elseif(NOT OGS_BUILD_TFEL)
        find_program(MFRONT mfront)
    endif()
    if(NOT MFRONT)
        BuildExternalProject(
            TFEL ${_tfel_source}
            CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                       "-DBUILD_SHARED_LIBS=OFF"
                       "-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
                       "-Denable-testing=OFF"
        )
        message(
            STATUS
                "ExternalProject_Add(): added package TFEL@rliv-${ogs.minimum_version.tfel-rliv}"
        )
        set(TFELHOME ${PROJECT_BINARY_DIR}/_ext/TFEL CACHE PATH "" FORCE)
        set(_EXT_LIBS ${_EXT_LIBS} TFEL CACHE INTERNAL "")
    endif()
endif()

if(OGS_USE_PETSC)
    option(OGS_BUILD_PETSC
           "Build PETSc locally. Needs to be set with a clean cache!" OFF
    )
    # Force CMake to accept a given PETSc configuration in case the failure of
    # MPI tests. This may cause the compilation broken.
    option(FORCE_PETSC_EXECUTABLE_RUNS
           "Force CMake to accept a given PETSc configuration" ON
    )
    if(FORCE_PETSC_EXECUTABLE_RUNS)
        set(PETSC_EXECUTABLE_RUNS YES)
    endif()

    set(_petsc_source GIT_REPOSITORY https://gitlab.com/petsc/petsc.git GIT_TAG
                      v${ogs.minimum_version.petsc}
    )
    set(_petsc_source_file
        ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/petsc-v${ogs.minimum_version.petsc}.zip
    )
    if(EXISTS ${_petsc_source_file})
        set(_petsc_source URL ${_petsc_source_file})
    elseif(NOT (OGS_PETSC_CONFIG_OPTIONS OR OGS_BUILD_PETSC))
        find_package(PETSc ${ogs.minimum_version.petsc})
    endif()

    if(NOT PETSC_FOUND)
        set(_configure_opts "")
        if(NOT "--download-fc=1" IN_LIST OGS_PETSC_CONFIG_OPTIONS)
            list(APPEND _configure_opts --with-fc=0)
        endif()
        if(ENV{CC})
            list(APPEND _configure_opts --with-cc=$ENV{CC})
        endif()
        if(ENV{CXX})
            list(APPEND _configure_opts --with-cxx=$ENV{CXX})
        endif()

        unset(ENV{PETSC_DIR})
        BuildExternalProject(
            PETSc ${_petsc_source}
            LOG_OUTPUT_ON_FAILURE ON
            CONFIGURE_COMMAND
                ./configure
                --download-f2cblaslapack=1
                --prefix=<INSTALL_DIR>
                --download-hypre
                --with-debugging=$<CONFIG:Debug>
                ${_configure_opts}
                ${OGS_PETSC_CONFIG_OPTIONS}
            BUILD_IN_SOURCE ON
            BUILD_COMMAND make -j all
            INSTALL_COMMAND make -j install
        )
        message(
            STATUS
                "ExternalProject_Add(): added package PETSc@${ogs.minimum_version.petsc}"
        )
        set(_EXT_LIBS ${_EXT_LIBS} PETSc CACHE INTERNAL "")
        BuildExternalProject_find_package(PETSc)
    endif()

    add_library(petsc SHARED IMPORTED)
    target_include_directories(petsc INTERFACE ${PETSC_INCLUDES})
    set_target_properties(petsc PROPERTIES IMPORTED_LOCATION ${PETSC_LIBRARIES})
    target_compile_definitions(petsc INTERFACE USE_PETSC)
endif()

if(OGS_USE_LIS)
    option(OGS_BUILD_LIS
           "Build LIS locally. Needs to be set with a clean cache!" OFF
    )
    set(_lis_source GIT_REPOSITORY https://github.com/anishida/lis.git GIT_TAG
                    ${ogs.minimum_version.lis}
    )
    set(_lis_source_file
        ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/lis-${ogs.minimum_version.lis}.zip
    )
    if(EXISTS ${_lis_source_file})
        set(_lis_source URL ${_lis_source_file})
    elseif(NOT OGS_BUILD_LIS)
        find_package(LIS)
    endif()
    if(NOT LIS_FOUND)
        if(BUILD_SHARED_LIBS)
            set(_lis_config_args --enable-shared)
        endif()
        BuildExternalProject(
            LIS ${_lis_source}
            CONFIGURE_COMMAND ./configure --enable-omp --prefix=<INSTALL_DIR>
                              ${_lis_config_args}
            BUILD_IN_SOURCE ON
            BUILD_COMMAND make -j
            INSTALL_COMMAND make -j install
        )
        message(
            STATUS
                "ExternalProject_Add(): added package LIS@${ogs.minimum_version.lis}"
        )
        set(ENV{LIS_ROOT_DIR} ${PROJECT_BINARY_DIR}/_ext/LIS)
        find_package(LIS REQUIRED)
    endif()
endif()

# ZLIB
option(OGS_BUILD_ZLIB "Build ZLIB locally. Needs to be set with a clean cache!"
       OFF
)
set(_zlib_source GIT_REPOSITORY https://github.com/madler/zlib.git GIT_TAG
                 v${ogs.tested_version.zlib}
)
set(_zlib_source_file
    ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/zlib-${ogs.tested_version.zlib}.zip
)
if(EXISTS ${_zlib_source_file})
    set(_zlib_source URL ${_zlib_source_file})
elseif(NOT OGS_BUILD_ZLIB)
    find_package(ZLIB)
endif()
if(NOT ZLIB_FOUND)
    BuildExternalProject(
        ZLIB ${_zlib_source}
        CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    )
    message(
        STATUS
            "ExternalProject_Add(): added package ZLIB@${ogs.tested_version.zlib}"
    )
    if(WIN32)
        # requires CMake 3.24 to be effective:
        set(ZLIB_USE_STATIC_LIBS "ON")
        set(ZLIB_ROOT ${PROJECT_BINARY_DIR}/_ext/ZLIB)
    endif()
    set(_EXT_LIBS ${_EXT_LIBS} ZLIB CACHE INTERNAL "")
    BuildExternalProject_find_package(ZLIB)
endif()

# HDF5
option(OGS_BUILD_HDF5 "Build HDF5 locally. Needs to be set with a clean cache!"
       OFF
)
set(_hdf5_options
    "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    "-DHDF5_GENERATE_HEADERS=OFF"
    "-DHDF5_BUILD_TOOLS=OFF"
    "-DHDF5_BUILD_EXAMPLES=OFF"
    "-DHDF5_BUILD_HL_LIB=${OGS_USE_PETSC}" # HL lib needed by MPI-enabled VTK
    "-DHDF5_BUILD_FORTRAN=OFF"
    "-DHDF5_BUILD_CPP_LIB=OFF"
    "-DHDF5_BUILD_JAVA=OFF"
    "-DBUILD_TESTING=OFF"
    "-DHDF5_ENABLE_Z_LIB_SUPPORT=ON"
)
if("${ZLIB_INCLUDE_DIRS}" MATCHES "${PROJECT_BINARY_DIR}/_ext/ZLIB")
    list(APPEND _hdf5_options "-DZLIB_ROOT=${PROJECT_BINARY_DIR}/_ext/ZLIB")
    if(WIN32)
        list(APPEND _hdf5_options "-DZLIB_USE_STATIC_LIBS=ON")
    endif()
endif()
if(OGS_USE_MPI)
    set(HDF5_PREFER_PARALLEL ON)
    list(APPEND _hdf5_options "-DHDF5_ENABLE_PARALLEL=ON")
endif()
if(WIN32 OR HDF5_USE_STATIC_LIBRARIES)
    set(HDF5_USE_STATIC_LIBRARIES ON)
    list(APPEND _hdf5_options "-DBUILD_SHARED_LIBS=OFF")
endif()
# Building from source requires newer hdf version
string(REPLACE "." "_" HDF5_TAG ${ogs.tested_version.hdf5})

set(_hdf5_source GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git GIT_TAG
                 hdf5-${HDF5_TAG}
)
set(_hdf5_source_file
    ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/hdf5-${ogs.tested_version.hdf5}.zip
)
if(EXISTS ${_hdf5_source_file})
    set(_hdf5_source URL ${_hdf5_source_file})
elseif(NOT OGS_BUILD_HDF5)
    find_package(HDF5 ${ogs.minimum_version.hdf5})
endif()
if(NOT HDF5_FOUND)
    BuildExternalProject(HDF5 ${_hdf5_source} CMAKE_ARGS ${_hdf5_options})
    message(
        STATUS
            "ExternalProject_Add(): added package HDF5@${ogs.tested_version.hdf5}"
    )
    set(_EXT_LIBS ${_EXT_LIBS} HDF5 CACHE INTERNAL "")
    BuildExternalProject_find_package(HDF5)
endif()

# append RPATHs
foreach(lib ${_EXT_LIBS})
    list(APPEND CMAKE_BUILD_RPATH ${PROJECT_BINARY_DIR}/_ext/${lib}/lib)
endforeach()
