get_source_files(SOURCES)
append_source_files(SOURCES IO)

append_source_files(SOURCES IO/XmlIO)

if(OGS_BUILD_GUI)
    append_source_files(SOURCES IO/XmlIO/Qt)
endif()

list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp" BuildInfo.h)

# Create the library
add_library(BaseLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS BaseLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

include(GenerateExportHeader)
generate_export_header(BaseLib)
target_include_directories(BaseLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(BaseLib PUBLIC logog)

if(MSVC)
    target_link_libraries(BaseLib PUBLIC WinMM) # needed for timeGetTime
endif()

if(OGS_BUILD_GUI)
    target_link_libraries(BaseLib PUBLIC Qt5::Xml Qt5::XmlPatterns)
endif()

if(OGS_USE_PCH)
    cotire(BaseLib)
endif()
