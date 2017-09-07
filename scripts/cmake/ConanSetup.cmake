if(NOT OGS_USE_CONAN)
    return()
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
    VTK/[>=7.1]@bilke/stable
    CACHE INTERNAL ""
)

set(CONAN_OPTIONS
    Boost:header_only=True
    Qt:xmlpatterns=True
    CACHE INTERNAL ""
)

if(LINUX AND BUILD_SHARED_LIBS)
    set(CONAN_OPTIONS ${CONAN_OPTIONS} VTK:fPIC=True)
endif()

if(OGS_BUILD_GUI)
    set(CONAN_REQUIRES ${CONAN_REQUIRES}
        Shapelib/1.3.0@bilke/stable
        libgeotiff/1.4.2@bilke/stable
        Qt/5.6.2@bilke/testing
    )
endif()

# Add Conan ogs remote
find_program(CONAN_CMD conan)
if(NOT CONAN_CMD)
    message(FATAL_ERROR "Conan executable not found!")
endif()
execute_process(COMMAND ${CONAN_CMD} remote list OUTPUT_VARIABLE CONAN_REMOTES)

# Add ogs remote
if("${CONAN_REMOTES}" MATCHES "ogs: https://ogs.jfrog.io/ogs/api/conan/conan")
    # Make sure ogs repo is first
    execute_process(COMMAND ${CONAN_CMD} remote update -i 0 ogs https://ogs.jfrog.io/ogs/api/conan/conan)
else()
    # Add ogs repo as first
    message(STATUS "Conan adding ogs remote repositoy (https://api.bintray.com/conan/ogs/conan)")
    execute_process(COMMAND ${CONAN_CMD} remote add -i 0 ogs https://ogs.jfrog.io/ogs/api/conan/conan)
endif()

# Add conan-community remote
if("${CONAN_REMOTES}" MATCHES "conan-community: https://api.bintray.com/conan/conan-community/conan")
    execute_process(COMMAND ${CONAN_CMD} remote update -i 2 conan-community https://api.bintray.com/conan/conan-community/conan)
else()
    message(STATUS "Conan adding community remote repositoy (https://api.bintray.com/conan/conan-community/conan)")
    execute_process(COMMAND ${CONAN_CMD} remote add -i 2 conan-community https://api.bintray.com/conan/conan-community/conan)
endif()

# Remove libraries from Conan which are set to "System"
message(STATUS "Third-party libraries:")
foreach(LIB ${OGS_LIBS})
    message("  - OGS_LIB_${LIB} = ${OGS_LIB_${LIB}}")
    if("${OGS_LIB_${LIB}}" STREQUAL System)
        list(FILTER CONAN_REQUIRES EXCLUDE REGEX ${LIB})
    endif()
endforeach()

conan_cmake_run(REQUIRES ${CONAN_REQUIRES}
                OPTIONS ${CONAN_OPTIONS}
                BASIC_SETUP
                UPDATE
                BUILD missing
)
