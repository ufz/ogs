if(NOT OGS_USE_CONAN)
    return()
endif()
string(TOLOWER ${OGS_USE_CONAN} OGS_USE_CONAN_lower)
if(OGS_USE_CONAN_lower STREQUAL "auto" AND POETRY)
    execute_process(COMMAND ${_CMD_COMMAND} poetry add conan==${ogs.minimum_version.conan}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    find_program(CONAN_CMD conan HINTS ${LOCAL_VIRTUALENV_BIN_DIRS}
        REQUIRED NO_DEFAULT_PATH
    )
else()
    find_program(CONAN_CMD conan)
endif()
if(NOT CONAN_CMD
   AND (OGS_USE_PETSC
        OR OGS_USE_LIS
        OR OGS_USE_MFRONT
        OR OGS_BUILD_GUI
       )
)
    message(WARNING "conan executable not found. Specify CMake option "
        "OGS_USE_CONAN=auto for automatic installation in the build directory "
        "OR install it system-wide (https://www.opengeosys.org/docs/devguide/"
        "getting-started/prerequisites/#step-install-conan-package-manager) "
        "OR disable this warning with OGS_USE_CONAN=OFF.")
    return()
endif()


if(CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(FATAL_ERROR "Multi-config generators are not yet supported when "
        "using Conan. Specify CMAKE_BUILD_TYPE, e.g. via cmd line: "
        "cmake . -DCMAKE_BUILD_TYPE=Release!")
endif()

# Treat Conan includes as system includes to suppress warnings
set(CONAN_SYSTEM_INCLUDES ON)

include(${PROJECT_SOURCE_DIR}/scripts/cmake/conan/conan.cmake)

if(OGS_USE_PETSC)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} petsc/${ogs.minimum_version.petsc}@bilke/testing)
    if(OGS_CONAN_USE_SYSTEM_OPENMPI)
        set(CONAN_OPTIONS ${CONAN_OPTIONS} petsc:skip_install_openmpi=True)
    endif()
endif()

if(OGS_USE_LIS)
    list(APPEND CONAN_OPTIONS lis:with_omp=True)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} lis/1.7.37@bilke/stable)
endif()

if(OGS_USE_MFRONT)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} tfel/3.4@bilke/testing)
endif()

if(OGS_BUILD_GUI)
    set(QT_VERSION ${ogs.minimum_version.qt})
    if(UNIX)
        set(QT_VERSION ${ogs.tested_version.qt})
    endif()
    set(CONAN_REQUIRES ${CONAN_REQUIRES}
        # libgeotiff/1.4.2@bilke/stable # TODO
        # Overrides for dependency mismatches
        bzip2/1.0.8
        zlib/1.2.11
        qt/${QT_VERSION}@bincrafters/stable
    )
    set(CONAN_OPTIONS ${CONAN_OPTIONS}
        qt:qtxmlpatterns=True
        qt:openssl=False
        qt:with_libalsa=False
        qt:with_libjpeg=False
        #qt:with_libpng=False
        qt:with_mysql=False
        qt:with_odbc=False
        qt:with_openal=False
        qt:with_pq=False
        qt:with_sdl2=False
        qt:with_sqlite3=False
    )
    if(MSVC)
        set(CONAN_OPTIONS ${CONAN_OPTIONS} qt:with_harfbuzz=False)
    endif()
    if(UNIX AND NOT APPLE)
        list(APPEND CONAN_OPTIONS qt:qtx11extras=True)
    endif()
endif()

if(OGS_USE_NETCDF)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} netcdf-cxx/4.3.1-1@bilke/testing)
endif()

if(NOT DEFINED CONAN_REQUIRES)
    return()
endif()

conan_check(VERSION ${ogs.minimum_version.conan})
conan_config_install(ITEM ${PROJECT_SOURCE_DIR}/scripts/cmake/conan/config)

message(STATUS "Third-party libraries:")
foreach(LIB ${OGS_LIBS})
    message(STATUS "  - OGS_LIB_${LIB} = ${OGS_LIB_${LIB}}")
    if("${OGS_LIB_${LIB}}" STREQUAL System)
        list(FILTER CONAN_REQUIRES EXCLUDE REGEX ${LIB})
    endif()
endforeach()

set(CONAN_IMPORTS "")
if(APPLE)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "lib, *.dylib* -> ./lib")
endif()
if(MSVC)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "bin, *.dll* -> ./bin")
endif()
if(UNIX AND NOT APPLE)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "lib, *.so* -> ./lib")
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "plugins/platforms, *.so* -> ./bin/platforms")
endif()

file(TIMESTAMP ${PROJECT_BINARY_DIR}/conan_install_timestamp.txt file_timestamp "%Y.%m.%d")
string(TIMESTAMP timestamp "%Y.%m.%d")

# Run conan install update only once a day
if("${file_timestamp}" VERSION_LESS ${timestamp} OR DEFINED ENV{CI})
    file(WRITE ${PROJECT_BINARY_DIR}/conan_install_timestamp.txt "${timestamp}\n")
    set(CONAN_UPDATE UPDATE)
    set(CONAN_COMMAND ${CONAN_CMD} CACHE INTERNAL "") # Speed up conan_add_remote
    conan_add_remote(NAME ogs INDEX 0
        URL https://ogs.jfrog.io/ogs/api/conan/conan)
    conan_add_remote(NAME bincrafters INDEX 1
        URL https://bincrafters.jfrog.io/artifactory/api/conan/public-conan)
    execute_process(COMMAND ${CONAN_COMMAND} config set general.revisions_enabled=1)
else()
    message(STATUS "Conan: Skipping update step.")
endif()

if(DEFINED OGS_CONAN_BUILD_TYPE)
    set(CONAN_BUILD_TYPE ${OGS_CONAN_BUILD_TYPE})
else()
    set(CONAN_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

if(MSVC)
    set(CC_CACHE $ENV{CC})
    set(CXX_CACHE $ENV{CXX})
    unset(ENV{CC}) # Disable clcache, e.g. for building qt
    unset(ENV{CXX})
endif()
# speed up conan_cmake_run
set(ARGUMENTS_CONAN_COMMAND ${CONAN_CMD} CACHE INTERNAL "")
conan_cmake_run(
    BASIC_SETUP
    ${CONAN_UPDATE}
    KEEP_RPATHS
    REQUIRES ${CONAN_REQUIRES}
    OPTIONS ${CONAN_OPTIONS}
    BUILD ${OGS_CONAN_BUILD}
    IMPORTS ${CONAN_IMPORTS}
    GENERATORS virtualrunenv
    BUILD_TYPE ${CONAN_BUILD_TYPE}
)
if(MSVC)
    set(ENV{CC} ${CC_CACHE}) # Restore vars
    set(ENV{CXX} ${CXX_CACHE})
endif()

if(NOT ${OGS_CONAN_BUILD} MATCHES "never|always|missing|outdated")
    message(STATUS "Warning: Resetting CMake variable OGS_CONAN_BUILD to its default value of 'missing'")
    set(OGS_CONAN_BUILD "missing" CACHE INTERNAL "")
endif()

if(OGS_USE_PETSC)
    set(PETSC_DIR ${CONAN_PETSC_ROOT} CACHE INTERNAL "")
endif()
