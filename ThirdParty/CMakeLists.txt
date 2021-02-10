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

if(OGS_BUILD_SWMM)
    add_subdirectory(SWMMInterface)
    list(APPEND DISABLE_WARNINGS_TARGETS SWMM SwmmInterface)
endif()

if(OGS_USE_XDMF)
    find_package(ZLIB REQUIRED) # ZLIB is a HDF5 dependency
    set(XDMF_LIBNAME OgsXdmf CACHE INTERNAL "")
    add_subdirectory(xdmf)
    target_include_directories(${XDMF_LIBNAME} PUBLIC
        ${PROJECT_SOURCE_DIR}/ThirdParty/xdmf
        ${PROJECT_BINARY_DIR}/ThirdParty/xdmf
    )
    if(OGS_USE_CONAN AND UNIX AND APPLE)
        find_package(Iconv REQUIRED)
    endif()

    if(MSVC AND OGS_USE_CONAN)
        # Hack: Conan HDF5 not found on Windows
        target_link_libraries(${XDMF_LIBNAME} ${CONAN_LIBS})
    else()
        target_link_libraries(${XDMF_LIBNAME} Boost::boost ${Iconv_LIBRARIES} ZLIB::ZLIB)
    endif()
    target_include_directories(${XDMF_LIBNAME}Core PUBLIC
        ${PROJECT_SOURCE_DIR}/ThirdParty/xdmf/core
        ${PROJECT_BINARY_DIR}/ThirdParty/xdmf/core
    )
    find_package(LibXml2 REQUIRED) # LibXml2 is a XdmfCore dependency
    target_link_libraries(${XDMF_LIBNAME}Core PUBLIC LibXml2::LibXml2)

    set_target_properties(${XDMF_LIBNAME} ${XDMF_LIBNAME}Core PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
    )
    if(BUILD_SHARED_LIBS)
        install(TARGETS ${XDMF_LIBNAME} ${XDMF_LIBNAME}Core LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
    list(APPEND DISABLE_WARNINGS_TARGETS ${XDMF_LIBNAME} ${XDMF_LIBNAME}Core)
endif()

foreach(TARGET ${DISABLE_WARNINGS_TARGETS})
    target_compile_options(${TARGET} PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-w>
        $<$<CXX_COMPILER_ID:MSVC>:/W0>)
endforeach()
