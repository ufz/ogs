add_subdirectory(spdlog)
set(DISABLE_WARNINGS_TARGETS spdlog)

if(BUILD_TESTING)
    set(gtest_force_shared_crt ON CACHE INTERNAL "") # Use dynamic MSVC runtime
    set(INSTALL_GTEST OFF CACHE INTERNAL "")
    add_subdirectory(googletest)
endif()

# This is a workaround for Travis builds.
option(OGS_USE_OPTIONAL_SUBMODULES "Option for enabling optional submodules" OFF)

# iphreeqc
if (OGS_BUILD_PROCESS_ComponentTransport
    OR OGS_BUILD_PROCESS_RichardsComponentTransport)
    add_subdirectory(iphreeqc)
    list(APPEND DISABLE_WARNINGS_TARGETS iphreeqc)
endif()

# MFrontGenericInterfaceSupport
if (OGS_USE_MFRONT)
    set(enable-doxygen-doc OFF CACHE INTERNAL "")
    set(enable-fortran-bindings OFF CACHE INTERNAL "")
    set(CMAKE_CXX_STANDARD 11)
    add_subdirectory(MGIS)
    set(CMAKE_CXX_STANDARD 17)
    list(APPEND DISABLE_WARNINGS_TARGETS MFrontGenericInterface)
endif()

# VtkFbxConverter
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VtkFbxConverter/CMakeLists.txt
    AND OGS_USE_OPTIONAL_SUBMODULES)
    add_subdirectory(VtkFbxConverter)
endif()

# tetgen
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tetgen/CMakeLists.txt)
    add_subdirectory(tetgen)
    list(APPEND DISABLE_WARNINGS_TARGETS tet tetgen)
    install(PROGRAMS $<TARGET_FILE:tetgen> DESTINATION bin COMPONENT ogs_extras)
endif()

if(OGS_BUILD_SWMM)
    add_subdirectory(SWMMInterface)
    list(APPEND DISABLE_WARNINGS_TARGETS SWMM SwmmInterface)
endif()

if(OGS_USE_PYTHON)
    # pybind11 uses old CMake find functionality, pass variables to use
    # the same Python installation.
    set(PYTHON_INCLUDE_DIR ${Python3_INCLUDE_DIRS} CACHE STRING "" FORCE)
    set(PYTHON_LIBRARIES ${Python3_LIBRARIES} CACHE STRING "" FORCE)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE} CACHE FILEPATH "" FORCE)
    set(PYBIND11_PYTHON_VERSION ${Python3_VERSION})

    add_subdirectory(pybind11)

    function(check_python_version_compatibility)
        if(NOT EXISTS ${VTK_DIR}/VTKConfig.cmake)
            return()
        endif()

        include(${VTK_DIR}/VTKConfig.cmake)
        if(NOT ${VTK_WRAP_PYTHON})
            return()
        endif()

        if(NOT EXISTS ${VTK_MODULES_DIR}/vtkPython.cmake)
            return()
        endif()

        include(${VTK_MODULES_DIR}/vtkPython.cmake)
        if (NOT "${vtkPython_LIBRARIES}" STREQUAL "${Python3_LIBRARIES}")
            message(SEND_ERROR "Mismatch between VTK's and OpenGeoSys' Python "
                "libraries: ${vtkPython_LIBRARIES} vs. ${Python3_LIBRARIES}. "
                "This will lead to compilation or linking errors. "
                "You can fix this error by using the same Python version for "
                "OpenGeoSys as VTK is built with.")
        endif()
    endfunction()

    check_python_version_compatibility()

endif()

if(OGS_BUILD_UTILS)
    include(${PROJECT_SOURCE_DIR}/scripts/cmake/MetisSetup.cmake)
    list(APPEND DISABLE_WARNINGS_TARGETS ogs_metis mpmetis)
endif()

foreach(TARGET ${DISABLE_WARNINGS_TARGETS})
    target_compile_options(${TARGET} PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-w>
        $<$<CXX_COMPILER_ID:MSVC>:/W0>)
endforeach()

find_package(ZLIB REQUIRED) # ZLIB is a HDF5 dependency
add_subdirectory(xdmf)
target_include_directories(Xdmf PUBLIC
    ${PROJECT_SOURCE_DIR}/ThirdParty/xdmf
    ${PROJECT_BINARY_DIR}/ThirdParty/xdmf
)
if(OGS_USE_CONAN AND UNIX AND APPLE)
    find_package(Iconv REQUIRED)
endif()
target_link_libraries(Xdmf Boost::boost ${Iconv_LIBRARIES} ZLIB::ZLIB)
target_include_directories(XdmfCore PUBLIC
    ${PROJECT_SOURCE_DIR}/ThirdParty/xdmf/core
    ${PROJECT_BINARY_DIR}/ThirdParty/xdmf/core
)
find_package(LibXml2 REQUIRED) # LibXml2 is a XdmfCore dependency
target_link_libraries(XdmfCore PUBLIC LibXml2::LibXml2)
