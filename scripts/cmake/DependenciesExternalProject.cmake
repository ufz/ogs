# because of ogs.minimum_version variables:
#
# cmake-lint: disable=C0103

# Build dependencies via ExternalProject_Add() at configure time in
# ${PROJECT_BINARY_DIR}/_ext or in $CPM_SOURCE_CACHE/_ext
include(BuildExternalProject)

set(OGS_EXTERNAL_DEPENDENCIES_CACHE ""
    CACHE PATH "Directory containing source archives of external dependencies."
)

set(_defaultCMakeArgs "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
if(CCACHE_EXECUTABLE)
    list(APPEND _defaultCMakeArgs "-DCMAKE_C_COMPILER_LAUNCHER=${CCACHE_EXECUTABLE}"
                                  "-DCMAKE_CXX_COMPILER_LAUNCHER=${CCACHE_EXECUTABLE}"
    )
endif()

if(MSVC)
    find_program(NINJA_CMD ninja)
    if(NINJA_CMD)
        # When ninja is available use it to speed up builds
        set(_cmake_generator CMAKE_GENERATOR Ninja)
        message(STATUS "Ninja generator will be used for external projects.")
    endif()
    # MKL env setup may override compiler
    list(APPEND _defaultCMakeArgs "-DCMAKE_CXX_COMPILER=cl.exe"
                                  "-DCMAKE_C_COMPILER=cl.exe"
    )
endif()

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
        if(CONDA_BUILD)
            set(TFELHOME $ENV{CONDA_PREFIX})
            set(CMAKE_REQUIRE_FIND_PACKAGE_TFEL TRUE)
        endif()
        find_package(TFEL)

        if(MFRONT AND APPLE)
            # TODO: check for version
            # ~~~
            # ➜ mfront --version
            # tfel
            # Version : 4.0.0
            # ~~~
            string(REPLACE "mfront" "" _mfront_bin_dir ${MFRONT})
            set(TFELHOME ${_mfront_bin_dir}/..)
        endif()
    endif()
    if(NOT MFRONT OR NOT TFEL_FOUND)
        set(_py_version_major_minor
            "${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}"
        )
        set(_py_boost_comp
            "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}"
        )
        find_package(Boost ${ogs.minimum_version.boost} COMPONENTS ${_py_boost_comp})
        if(Boost_${_py_boost_comp}_FOUND)
            set(_tfel_cmake_args
                "-DPython_ADDITIONAL_VERSIONS=${_py_version_major_minor}"
                "-Denable-python-bindings=ON"
                "-Denable-numpy-support=OFF"
            )
        else()
            # Cleanup variables from previous find_package()-call
            unset(Boost_INCLUDE_DIR)
            unset(Boost_INCLUDE_DIRS)
            message(
                STATUS
                    "TFEL Python bindings disabled as Boosts Python library was not found."
            )
        endif()

        # Only one flag supported, prefer ASAN
        if(ENABLE_ASAN)
            set(_sanitize_flag -fsanitize=address)
        endif()
        if(ENABLE_UBSAN AND NOT DEFINED _sanitize_flag)
            set(_sanitize_flag -fsanitize=undefined)
        elseif(ENABLE_UBSAN AND DEFINED _sanitize_flag)
            message(STATUS "MFront: ASAN enabled only! UBSAN is off.")
        endif()
        if(DEFINED _sanitize_flag)
            foreach(var CXX EXE_LINKER SHARED_LINKER MODULE_LINKER)
                list(APPEND _tfel_cmake_args
                     "-DCMAKE_${var}_FLAGS_INIT=${_sanitize_flag}"
                )
            endforeach()
        endif()

        BuildExternalProject(
            TFEL ${_tfel_source}
            CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                       "-DBUILD_SHARED_LIBS=OFF"
                       "-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
                       "-Denable-testing=OFF"
                       ${_defaultCMakeArgs}
                       "${_tfel_cmake_args}"
                       SKIP_FIND
        )
        message(
            STATUS
                "ExternalProject_Add(): added package TFEL@rliv-${ogs.minimum_version.tfel-rliv}"
        )
        if(Boost_${_py_boost_comp}_FOUND)
            set(TFEL_WITH_PYTHON
                "${build_dir_TFEL}/lib/python${_py_version_major_minor}/site-packages"
                CACHE PATH ""
            )
        endif()
        set(TFELHOME ${build_dir_TFEL} CACHE PATH "" FORCE)
    endif()
    if(TFEL_WITH_PYTHON)
        message(STATUS "TFEL build with Python bindings. To use them:\n "
                       "  export PYTHONPATH=${TFEL_WITH_PYTHON}:$PYTHONPATH"
        )
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

    # apple clang 15 requires newer petsc, see
    # https://www.mail-archive.com/petsc-users@mcs.anl.gov/msg46980.html
    if(APPLE AND ${CMAKE_CXX_COMPILER_VERSION} GREATER_EQUAL 15)
        if(3.20.5 VERSION_GREATER ${ogs.minimum_version.petsc})
            set(ogs.minimum_version.petsc 3.20.5)
        endif()
    endif()

    # On Python 3.13 this fix is needed:
    # https://gitlab.com/petsc/petsc/-/merge_requests/7790
    if(Python_VERSION GREATER_EQUAL 3.13)
        if(3.22.4 VERSION_GREATER ${ogs.minimum_version.petsc})
            set(ogs.minimum_version.petsc 3.22.4)
        endif()
    endif()

    set(_petsc_source GIT_REPOSITORY https://gitlab.com/petsc/petsc.git GIT_TAG
                      v${ogs.minimum_version.petsc}
    )
    set(_petsc_source_file
        ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/petsc-v${ogs.minimum_version.petsc}.zip
    )
    if(DEFINED ENV{OGS_PETSC_CONFIG_OPTIONS} AND "${OGS_PETSC_CONFIG_OPTIONS}"
                                                 STREQUAL ""
    )
        set(OGS_PETSC_CONFIG_OPTIONS "$ENV{OGS_PETSC_CONFIG_OPTIONS}")
    endif()
    if(EXISTS ${_petsc_source_file})
        set(_petsc_source URL ${_petsc_source_file})
    elseif(NOT (OGS_PETSC_CONFIG_OPTIONS OR OGS_BUILD_PETSC))
        find_package(PkgConfig REQUIRED)
        pkg_search_module(PETSc IMPORTED_TARGET PETSc)
    endif()

    if(NOT TARGET PkgConfig::PETSc)
        set(_configure_opts "")
        if(NOT "--download-fc" IN_LIST OGS_PETSC_CONFIG_OPTIONS)
            list(APPEND _configure_opts --with-fc=0)
        endif()

        unset(ENV{PETSC_DIR})
        BuildExternalProject(
            PETSc ${_petsc_source}
            LOG_OUTPUT_ON_FAILURE ON
            CONFIGURE_COMMAND
                ./configure
                --prefix=<INSTALL_DIR>
                --download-f2cblaslapack=1
                --download-hypre=1
                --with-debugging=$<CONFIG:Debug>
                ${_configure_opts}
                ${OGS_PETSC_CONFIG_OPTIONS}
            BUILD_IN_SOURCE ON
            BUILD_COMMAND make -j$ENV{CMAKE_BUILD_PARALLEL_LEVEL} all
            INSTALL_COMMAND make -j$ENV{CMAKE_BUILD_PARALLEL_LEVEL} install
                            SKIP_FIND
        )
        message(
            STATUS
                "ExternalProject_Add(): added package PETSc@${ogs.minimum_version.petsc}"
        )
        set(_EXT_LIBS ${_EXT_LIBS} PETSc CACHE INTERNAL "")
        set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${build_dir_PETSc})
    endif()

    find_package(PkgConfig REQUIRED)
    pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
    target_compile_definitions(PkgConfig::PETSC INTERFACE USE_PETSC)
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
            BUILD_COMMAND make -j$ENV{CMAKE_BUILD_PARALLEL_LEVEL}
            INSTALL_COMMAND make -j$ENV{CMAKE_BUILD_PARALLEL_LEVEL} install
        )
        message(
            STATUS
                "ExternalProject_Add(): added package LIS@${ogs.minimum_version.lis}"
        )
        set(ENV{LIS_ROOT_DIR} ${build_dir_LIS})
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
        ZLIB ${_zlib_source} CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                                        "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
                                        ${_defaultCMakeArgs}
    )
    message(
        STATUS
            "ExternalProject_Add(): added package ZLIB@${ogs.tested_version.zlib}"
    )
    if(WIN32)
        if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(_zlib_debug_postfix "d")
        endif()
        # Remove zlib dll and corresponding lib, ZLIB_USE_STATIC_LIBS sometimes
        # does not work
        file(
            REMOVE
            ${build_dir_ZLIB}/${CMAKE_INSTALL_LIBDIR}/zlib${_zlib_debug_postfix}${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${build_dir_ZLIB}/${CMAKE_INSTALL_BINDIR}/zlib${_zlib_debug_postfix}${CMAKE_SHARED_LIBRARY_SUFFIX}
        )
        # requires CMake 3.24 to be effective:
        set(ZLIB_USE_STATIC_LIBS "ON")
        set(ZLIB_ROOT ${build_dir_ZLIB})
        # Force local zlib build, found netcdf-installed zlib sometimes
        set(ZLIB_LIBRARIES
            ${build_dir_ZLIB}/${CMAKE_INSTALL_LIBDIR}/zlibstatic${_zlib_debug_postfix}${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
        set(ZLIB_INCLUDE_DIRS ${build_dir_ZLIB})
    endif()
    set(_EXT_LIBS ${_EXT_LIBS} ZLIB CACHE INTERNAL "")
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
if("${ZLIB_INCLUDE_DIRS}" MATCHES "${build_dir_ZLIB}")
    list(APPEND _hdf5_options "-DZLIB_ROOT=${build_dir_ZLIB}")
    if(WIN32)
        list(APPEND _hdf5_options "-DZLIB_USE_STATIC_LIBS=ON")
    endif()
endif()
if(OGS_USE_MPI)
    list(APPEND _hdf5_options "-DHDF5_ENABLE_PARALLEL=ON")
endif()
if(WIN32 OR HDF5_USE_STATIC_LIBRARIES)
    set(HDF5_USE_STATIC_LIBRARIES ON)
    list(APPEND _hdf5_options "-DBUILD_SHARED_LIBS=OFF")
endif()

# With apple clang 15 there are errors when the cpm compiled hdf5 has a
# different version than the one bundled with vtk, see
# https://gitlab.kitware.com/vtk/vtk/-/issues/19232.
if(APPLE AND ${CMAKE_CXX_COMPILER_VERSION} GREATER_EQUAL 15)
    set(ogs.tested_version.hdf5 1.10.7)
endif()

# Building from source requires newer hdf version
string(REPLACE "." "_" HDF5_TAG ${ogs.tested_version.hdf5})

set(_hdf5_source GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git GIT_TAG
                 hdf5-${HDF5_TAG}
)
set(_hdf5_source_file
    ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/hdf5-${ogs.tested_version.hdf5}.zip
)

set(HDF5_PREFER_PARALLEL ON)

if(EXISTS ${_hdf5_source_file})
    set(_hdf5_source URL ${_hdf5_source_file})
elseif(NOT OGS_BUILD_HDF5)
    find_package(HDF5 ${ogs.minimum_version.hdf5} COMPONENTS C)
endif()
if(NOT _HDF5_FOUND AND NOT HDF5_FOUND)
    BuildExternalProject(
        HDF5 ${_hdf5_source} CMAKE_ARGS ${_hdf5_options} ${_defaultCMakeArgs}
                                        ${_cmake_generator} SKIP_FIND
    )
    message(
        STATUS
            "ExternalProject_Add(): added package HDF5@${ogs.tested_version.hdf5}"
    )
    set(_EXT_LIBS ${_EXT_LIBS} HDF5 CACHE INTERNAL "")
    set(_HDF5_FOUND ON CACHE INTERNAL "")
endif()
if(_HDF5_FOUND)
    # Use HDF-provided install method, see ext build directory /
    # share/USING_HDF5_CMake.txt
    if(HDF5_USE_STATIC_LIBRARIES)
        set(_hdf5_LIB_TYPE STATIC) # or SHARED
    else()
        set(_hdf5_LIB_TYPE SHARED) # or SHARED
    endif()
    string(TOLOWER ${_hdf5_LIB_TYPE} _hdf5_SEARCH_TYPE)
    set(HDF5_ROOT ${build_dir_HDF5})
    find_package(
        HDF5 REQUIRED NAMES hdf5
        COMPONENTS C ${_hdf5_SEARCH_TYPE} NO_CMAKE_PATH
                   NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH
                   NO_CMAKE_SYSTEM_PATH
    )
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
    set(HDF5_LIBRARIES ${HDF5_C_${_hdf5_LIB_TYPE}_LIBRARY})
    set(HDF5_C_LIBRARIES ${HDF5_LIBRARIES})
    if(WIN32)
        set(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${ZLIB_LIBRARIES})
    endif()
    message(STATUS "Using HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
endif()

# VTK
option(OGS_BUILD_VTK "Build VTK locally. Needs to be set with a clean cache!"
       OFF
)
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
            list(
                APPEND
                VTK_OPTIONS
                "-D${ogs.libraries.vtk.options_${option_index}.cmake_${cmake_index}}"
            )
        endforeach()
    endif()
endforeach()
list(REMOVE_DUPLICATES VTK_OPTIONS)

# Setting static libs for easier packaging.
list(APPEND VTK_OPTIONS "-DBUILD_SHARED_LIBS=OFF"
     "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
)
message(STATUS "VTK_OPTIONS: ${VTK_OPTIONS}")
if(OGS_USE_PETSC AND EXISTS ${build_dir_HDF5})
    # Use local hdf5 build
    list(APPEND VTK_OPTIONS "-DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON"
         "-DHDF5_ROOT=${build_dir_HDF5}"
    )
endif()

# Building from source requires newer hdf version
string(REPLACE "." "_" HDF5_TAG ${ogs.tested_version.hdf5})

set(_vtk_source GIT_REPOSITORY https://github.com/kitware/vtk.git GIT_TAG
                v${ogs.minimum_version.vtk}
)
set(_vtk_source_file
    ${OGS_EXTERNAL_DEPENDENCIES_CACHE}/vtk-v${ogs.minimum_version.vtk}.tar.gz
)
if(GUIX_BUILD)
    find_package(VTK COMPONENTS ${VTK_COMPONENTS})
elseif(EXISTS ${_vtk_source_file})
    set(_vtk_source URL ${_vtk_source_file})
elseif(NOT OGS_BUILD_VTK AND (NOT OGS_USE_MKL OR CONDA_BUILD))
    # Typically VTK also pulls in libgomp dependency when found on system
    unset(VTK_COMPONENTS)
    foreach(opt ${VTK_OPTIONS})
        if("${opt}" MATCHES "-DVTK_MODULE_ENABLE_VTK_(.*)=YES")
            list(APPEND VTK_COMPONENTS ${CMAKE_MATCH_1})
        endif()
    endforeach()
    message(STATUS "Searching VTK on system with components: ${VTK_COMPONENTS}")
    find_package(VTK ${ogs.minimum_version.vtk} COMPONENTS ${VTK_COMPONENTS})
endif()
if(NOT VTK_FOUND)
    file(
        DOWNLOAD
        https://gitlab.kitware.com/bilke/vtk/-/commit/b70e3e103cf711e080f23171201c7d030187146b.patch
        ${PROJECT_SOURCE_DIR}/scripts/cmake/vtk-win.patch
    )
    file(
        DOWNLOAD
        https://gitlab.kitware.com/bilke/vtk/-/commit/70b16fda87f82520fa29b48c6a62bafa405d8ee2.patch
        ${PROJECT_SOURCE_DIR}/scripts/cmake/vtk-mac.patch
    )
    if("${OGS_EXTERNAL_DEPENDENCIES_CACHE}" STREQUAL "")
        set(_vtk_patch PATCH_COMMAND git apply)
        if(WIN32)
            # Fixes https://gitlab.kitware.com/vtk/vtk/-/issues/19178
            list(APPEND _vtk_patch
                 "${PROJECT_SOURCE_DIR}/scripts/cmake/vtk-win.patch"
            )
            message(STATUS "Applying VTK Win patch")
        endif()
        if(APPLE)
            # Fixes https://stackoverflow.com/questions/9894961
            list(APPEND _vtk_patch
                 "${PROJECT_SOURCE_DIR}/scripts/cmake/vtk-mac.patch"
            )
            message(STATUS "Applying VTK Mac patch")
        endif()
        if(LINUX)
            # No patches on Linux
            unset(_vtk_patch)
        endif()
    endif()

    if(DEFINED Qt5_DIR)
        list(APPEND VTK_OPTIONS "-DQt5_DIR=${Qt5_DIR}")
    endif()

    BuildExternalProject(
        VTK ${_vtk_source} CMAKE_ARGS ${VTK_OPTIONS} ${_defaultCMakeArgs}
                                      ${_vtk_patch} ${_cmake_generator}
    )
    message(
        STATUS
            "ExternalProject_Add(): added package VTK@${ogs.minimum_version.vtk}"
    )
    set(VTK_VERSION "${ogs.minimum_version.vtk}" CACHE INTERNAL "")
    set(_EXT_LIBS ${_EXT_LIBS} VTK CACHE INTERNAL "")
endif()

# cmake-lint: disable=C0103

# append RPATHs
foreach(lib ${_EXT_LIBS})
    set(CMAKE_BUILD_RPATH ${CMAKE_BUILD_RPATH} ${build_dir_${lib}}/lib
                          ${build_dir_${lib}}/lib64
    )
    set(${lib}_SOURCE_DIR ${build_dir_${lib}}/src/${lib})
endforeach()
