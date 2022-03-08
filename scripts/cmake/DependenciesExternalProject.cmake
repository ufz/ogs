# Build dependencies via ExternalProject_Add() at configure time in
# ${PROJECT_BINARY_DIR}/_ext
include(BuildExternalProject)

set(OGS_EXTERNAL_DEPENDENCIES_CACHE ""
    CACHE PATH "Directory containing source archives of external dependencies."
)

if(OGS_INSTALL_EXTERNAL_DEPENDENCIES)
    set(_install_dir INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
endif()

if(OGS_USE_MFRONT)
    set(_tfel_source GIT_REPOSITORY https://github.com/thelfer/tfel.git GIT_TAG
                     rliv-${ogs.minimum_version.tfel-rliv}
    )
    set(_tfel_source_file
        ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/tfel-rliv-${ogs.minimum_version.tfel-rliv}.zip
    )
    if(EXISTS ${_tfel_source_file})
        set(_tfel_source URL ${_tfel_source_file})
    else()
        find_program(MFRONT mfront)
    endif()
    if(NOT MFRONT)
        BuildExternalProject(
            TFEL ${_tfel_source} ${_install_dir}
            CMAKE_ARGS -DCMAKE_INSTALL_RPATH=<INSTALL_DIR>/lib
                       -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE
        )
        message(
            STATUS
                "ExternalProject_Add(): added package TFEL@rliv-${ogs.minimum_version.tfel-rliv}"
        )
        if(OGS_INSTALL_EXTERNAL_DEPENDENCIES)
            set(TFELHOME ${CMAKE_INSTALL_PREFIX} CACHE PATH "" FORCE)
        else()
            set(TFELHOME ${PROJECT_BINARY_DIR}/_ext/TFEL CACHE PATH "" FORCE)
        endif()
    endif()
    list(APPEND CMAKE_INSTALL_RPATH ${TFELHOME}/lib)
endif()

if(OGS_USE_PETSC)
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
    elseif(NOT OGS_PETSC_CONFIG_OPTIONS)
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
                ./configure --download-f2cblaslapack=1 --prefix=<INSTALL_DIR>
                --with-debugging=$<CONFIG:Debug> ${_configure_opts}
                ${OGS_PETSC_CONFIG_OPTIONS}
            BUILD_IN_SOURCE ON
            BUILD_COMMAND make -j all
            INSTALL_COMMAND make -j install ${_install_dir}
        )
        message(
            STATUS
                "ExternalProject_Add(): added package PETSc@${ogs.minimum_version.petsc}"
        )
        if(OGS_INSTALL_EXTERNAL_DEPENDENCIES)
            set(PETSC_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "" FORCE)
        else()
            set(PETSC_DIR ${PROJECT_BINARY_DIR}/_ext/PETSc CACHE PATH "" FORCE)
        endif()
        find_package(PETSc ${ogs.minimum_version.petsc} REQUIRED)
    endif()

    add_library(petsc SHARED IMPORTED)
    target_include_directories(petsc INTERFACE ${PETSC_INCLUDES})
    set_target_properties(petsc PROPERTIES IMPORTED_LOCATION ${PETSC_LIBRARIES})
    target_compile_definitions(petsc INTERFACE USE_PETSC)

    if(EXISTS ${PETSC_DIR}/lib)
        message(STATUS "RPATH: Appending ${PETSC_DIR}/lib")
        list(APPEND CMAKE_INSTALL_RPATH ${PETSC_DIR}/lib)
    endif()
endif()

if(OGS_USE_LIS)
    set(_lis_source GIT_REPOSITORY https://github.com/anishida/lis.git GIT_TAG
                    ${ogs.minimum_version.lis}
    )
    set(_lis_source_file
        ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/lis-${ogs.minimum_version.lis}.zip
    )
    if(EXISTS ${_lis_source_file})
        set(_lis_source URL ${_lis_source_file})
    else()
        find_package(LIS)
    endif()
    if(NOT LIS_FOUND)
        BuildExternalProject(
            LIS ${_lis_source}
            CONFIGURE_COMMAND ./configure --enable-omp --prefix=<INSTALL_DIR>
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
