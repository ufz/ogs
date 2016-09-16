# Source files
GET_SOURCE_FILES(SOURCES)

APPEND_SOURCE_FILES(SOURCES IO)
APPEND_SOURCE_FILES(SOURCES IO/Legacy)

APPEND_SOURCE_FILES(SOURCES IO/XmlIO/Rapid)
APPEND_SOURCE_FILES(SOURCES IO/XmlIO/Boost)

if(QT4_FOUND)
    APPEND_SOURCE_FILES(SOURCES IO/XmlIO/Qt)
endif()

# Create the library
add_library(GeoLib ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/tetgen/predicates.cxx
)

target_link_libraries(GeoLib
    MathLib
)

if(QT4_FOUND)
    target_link_libraries(GeoLib Qt4::QtXml Qt4::QtXmlPatterns)
    if(WIN32 AND CMAKE_CROSSCOMPILING AND OPENSSL_FOUND)
        target_link_libraries(GeoLib Qt4::QtNetwork ${OPENSSL_LIBRARIES} ws2_32)
    endif()
endif()

if(TARGET Boost)
    add_dependencies(GeoLib Boost)
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/packaging/InstallXmlSchemaFiles.cmake)
InstallXmlSchemaFiles("IO/XmlIO/*.xsd")
