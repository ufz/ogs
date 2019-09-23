if(NOT OGS_USE_CONAN)
    return()
endif()
find_program(CONAN_CMD conan)
if(NOT CONAN_CMD)
    message(WARNING "conan executable not found. Consider installing Conan for "
        "automatic third-party library handling. https://www.opengeosys.org/doc"
        "s/devguide/getting-started/prerequisites/#step-install-conan-package-m"
        "anager OR disable this warning with OGS_USE_CONAN=OFF")
    return()
endif()


if(CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(FATAL_ERROR "Multi-config generators are not yet supported when "
        "using Conan. Specify CMAKE_BUILD_TYPE!")
endif()

# $ cat /etc/os-release | grep VERSION_ID
# VERSION_ID="17.10"
if(COMPILER_IS_GCC AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 6 AND UBUNTU_VERSION VERSION_EQUAL 16)

endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/conan/conan.cmake)

set(CONAN_REQUIRES
    boost/1.66.0@conan/stable
    eigen/3.3.4@conan/stable
    vtk/8.2.0@bilke/stable
    CACHE INTERNAL ""
)

set(CONAN_OPTIONS
    boost:header_only=True
    vtk:minimal=True
    vtk:ioxml=True
    CACHE INTERNAL ""
)

if((UNIX AND NOT APPLE) AND BUILD_SHARED_LIBS)
    set(CONAN_OPTIONS ${CONAN_OPTIONS} vtk:fPIC=True)
endif()

if(OGS_USE_MPI)
    set(CONAN_OPTIONS ${CONAN_OPTIONS} vtk:mpi_minimal=True)
endif()

if(OGS_USE_PETSC)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} petsc/3.8.3@bilke/testing)
    if(OGS_CONAN_USE_SYSTEM_OPENMPI)
        set(CONAN_OPTIONS ${CONAN_OPTIONS} petsc:skip_install_openmpi=True)
    endif()
endif()

if(OGS_USE_LIS)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} lis/1.7.9@bilke/stable)
endif()

if(OGS_USE_CVODE)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} cvode/2.8.2@bilke/stable)
endif()

if(OGS_USE_MFRONT)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} tfel/3.2.1@bilke/testing)
endif()

if(OGS_BUILD_GUI)
    set(CONAN_REQUIRES ${CONAN_REQUIRES}
        shapelib/1.3.0@bilke/stable
        libgeotiff/1.4.2@bilke/stable
        qt/5.12.4@bincrafters/stable
        # Overwrite VTK requirement to match Qt requirement
        bzip2/1.0.8@conan/stable
    )
    set(CONAN_OPTIONS ${CONAN_OPTIONS}
        vtk:minimal=False
        vtk:qt=True
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
endif()

if(OGS_USE_NETCDF)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} netcdf-cxx/4.3.1@bilke/testing)
endif()

conan_check(VERSION 1.3.0)

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
if("${file_timestamp}" VERSION_LESS ${timestamp} OR IS_CI)
    file(WRITE ${PROJECT_BINARY_DIR}/conan_install_timestamp.txt "${timestamp}\n")
    set(CONAN_UPDATE UPDATE)
    conan_add_remote(NAME ogs INDEX 0
        URL https://ogs.jfrog.io/ogs/api/conan/conan)
    conan_add_remote(NAME conan-community INDEX 1
        URL https://api.bintray.com/conan/conan-community/conan)
    conan_add_remote(NAME bincrafters INDEX 2
        URL https://api.bintray.com/conan/bincrafters/public-conan)
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
