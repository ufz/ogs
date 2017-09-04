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
    Boost/[>=1.56.0]@lasote/stable
    Eigen3/3.2.9@bilke/stable
    VTK/[>=7.1]@bilke/stable
    CACHE INTERNAL ""
)

set(CONAN_OPTIONS
    Boost:header_only=True
    Boost:without_iostreams=True
    Qt:xmlpatterns=True
    CACHE INTERNAL ""
)

if(OGS_BUILD_GUI)
    set(CONAN_REQUIRES ${CONAN_REQUIRES}
        Shapelib/1.3.0@bilke/stable
        libgeotiff/1.4.2@bilke/stable
        Qt/5.6.2@bilke/testing
    )
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
                BUILD missing
)
