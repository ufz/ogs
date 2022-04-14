# Source files
get_source_files(SOURCES)

append_source_files(SOURCES IO)

append_source_files(SOURCES IO/XmlIO/Boost)

if(OGS_BUILD_GUI)
    append_source_files(SOURCES IO/XmlIO/Qt)
endif()

# Create the library
ogs_add_library(GeoLib ${SOURCES})

target_link_libraries(
    GeoLib
    PUBLIC BaseLib Eigen3::Eigen MathLib $<$<BOOL:${OGS_BUILD_GUI}>:rapidxml>
           $<$<BOOL:${OGS_BUILD_GUI}>:Qt5::Xml>
           $<$<BOOL:${OGS_BUILD_GUI}>:Qt5::XmlPatterns>
    PRIVATE tet
)

foreach(xsd OpenGeoSysGLI.xsd OpenGeoSysSTN.xsd)
    # cmake-lint: disable=E1126
    file(COPY_FILE IO/XmlIO/${xsd} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${xsd}
         ONLY_IF_DIFFERENT
    )
    install(FILES IO/XmlIO/${xsd} DESTINATION bin)
endforeach()
