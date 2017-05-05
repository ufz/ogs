# Source files
GET_SOURCE_FILES(SOURCES)

APPEND_SOURCE_FILES(SOURCES IO)
APPEND_SOURCE_FILES(SOURCES IO/Legacy)

APPEND_SOURCE_FILES(SOURCES IO/XmlIO/Rapid)
APPEND_SOURCE_FILES(SOURCES IO/XmlIO/Boost)

if(Qt5XmlPatterns_FOUND)
    APPEND_SOURCE_FILES(SOURCES IO/XmlIO/Qt)
endif()

# Create the library
add_library(GeoLib ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/tetgen/predicates.cxx
)

target_link_libraries(GeoLib PUBLIC BaseLib MathLib logog)

if(Qt5XmlPatterns_FOUND)
    target_link_libraries(GeoLib PUBLIC Qt5::Xml Qt5::XmlPatterns)
    if(WIN32 AND CMAKE_CROSSCOMPILING AND OPENSSL_FOUND)
        target_link_libraries(GeoLib PUBLIC Qt5::Network ${OPENSSL_LIBRARIES} ws2_32)
    endif()
endif()

if(TARGET Boost)
    add_dependencies(GeoLib Boost)
endif()

if(OGS_USE_PCH)
    cotire(GeoLib)
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/packaging/InstallXmlSchemaFiles.cmake)
InstallXmlSchemaFiles("IO/XmlIO/*.xsd")
