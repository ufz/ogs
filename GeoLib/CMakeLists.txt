# Source files
get_source_files(SOURCES)

append_source_files(SOURCES IO)

append_source_files(SOURCES IO/XmlIO/Rapid)
append_source_files(SOURCES IO/XmlIO/Boost)

if(OGS_BUILD_GUI)
    append_source_files(SOURCES IO/XmlIO/Qt)
endif()

# Create the library
add_library(GeoLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS GeoLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(GeoLib PUBLIC BaseLib MathLib logog PRIVATE tet)

if(OGS_BUILD_GUI)
    target_link_libraries(GeoLib PUBLIC Qt5::Xml Qt5::XmlPatterns)
endif()

if(OGS_USE_PCH)
    cotire(GeoLib)
endif()

configure_file(IO/XmlIO/OpenGeoSysGLI.xsd
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/OpenGeoSysGLI.xsd COPYONLY)
configure_file(IO/XmlIO/OpenGeoSysSTN.xsd
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/OpenGeoSysSTN.xsd COPYONLY)
install(FILES IO/XmlIO/OpenGeoSysGLI.xsd IO/XmlIO/OpenGeoSysSTN.xsd
        DESTINATION bin
        COMPONENT ogs_cli)
