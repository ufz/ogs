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

if(DEFINED OGS_LIB_Boost)
    if(${OGS_LIB_Boost} STREQUAL "Default")
        cmake_minimum_required(VERSION 3.4) # Conan Boost package requires this
    endif()
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/conan/conan.cmake)

set(CONAN_REQUIRES
    Boost/1.64.0@conan/stable
    Eigen3/3.2.9@bilke/stable
    VTK/8.1.1@bilke/stable
    CACHE INTERNAL ""
)

set(CONAN_OPTIONS
    Boost:header_only=True
    Qt:qtxmlpatterns=True
    VTK:minimal=True
    VTK:ioxml=True
    CACHE INTERNAL ""
)

if((UNIX AND NOT APPLE) AND BUILD_SHARED_LIBS)
    set(CONAN_OPTIONS ${CONAN_OPTIONS} VTK:fPIC=True)
endif()

if(OGS_USE_MPI)
    set(CONAN_OPTIONS ${CONAN_OPTIONS} VTK:mpi_minimal=True)
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

if(OGS_BUILD_GUI)
    set(CONAN_REQUIRES ${CONAN_REQUIRES}
        Shapelib/1.3.0@bilke/stable
        libgeotiff/1.4.2@bilke/stable
        Qt/5.11.0@bilke/stable
    )
    set(CONAN_OPTIONS ${CONAN_OPTIONS}
        VTK:minimal=False
        VTK:qt=True
    )
endif()

conan_check(VERSION 1.0.0)
conan_add_remote(NAME ogs INDEX 0
    URL https://ogs.jfrog.io/ogs/api/conan/conan)
conan_add_remote(NAME conan-community INDEX 1
    URL https://api.bintray.com/conan/conan-community/conan)
conan_add_remote(NAME bincrafters INDEX 2
    URL https://api.bintray.com/conan/bincrafters/public-conan)

message(STATUS "Third-party libraries:")
foreach(LIB ${OGS_LIBS})
    message("  - OGS_LIB_${LIB} = ${OGS_LIB_${LIB}}")
    if("${OGS_LIB_${LIB}}" STREQUAL System)
        list(FILTER CONAN_REQUIRES EXCLUDE REGEX ${LIB})
    endif()
endforeach()

set(CONAN_IMPORTS "")
if(APPLE)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "lib, *.dylib* -> ./bin")
endif()
if(MSVC)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "bin, *.dll* -> ./bin")
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "plugins/platforms, *.dll* -> ./bin/platforms")
endif()
if(UNIX AND NOT APPLE)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "lib, *.so* -> ./lib")
endif()

file(TIMESTAMP ${PROJECT_BINARY_DIR}/conan_install_timestamp.txt file_timestamp "%Y.%m.%d")
string(TIMESTAMP timestamp "%Y.%m.%d")

# Run conan install update only once a day
if("${file_timestamp}" VERSION_LESS ${timestamp} OR IS_CI)
    file(WRITE ${PROJECT_BINARY_DIR}/conan_install_timestamp.txt "${timestamp}\n")
    set(CONAN_UPDATE UPDATE)
else()
    message(STATUS "Conan: Skipping update step.")
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
)

if(NOT ${OGS_CONAN_BUILD} MATCHES "never|always|missing")
    message(STATUS "Warning: Resetting CMake variable OGS_CONAN_BUILD to its default value of 'missing'")
    set(OGS_CONAN_BUILD "missing" CACHE INTERNAL "")
endif()

if(OGS_USE_PETSC)
    set(PETSC_DIR ${CONAN_PETSC_ROOT} CACHE INTERNAL "")
endif()
