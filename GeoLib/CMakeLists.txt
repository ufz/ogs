# Source files
GET_SOURCE_FILES(SOURCES_GeoLib)
set(SOURCES ${SOURCES_GeoLib})

GET_SOURCE_FILES(SOURCES_IO IO)
GET_SOURCE_FILES(SOURCES_IO_Legacy IO/Legacy)
set(SOURCES ${SOURCES} ${SOURCES_IO} ${SOURCES_IO_Legacy})

GET_SOURCE_FILES(SOURCES_IO_BOOST_XML IO/XmlIO/Boost)
GET_SOURCE_FILES(SOURCES_IO_RAPID_XML IO/XmlIO/Rapid)
set(SOURCES ${SOURCES} ${SOURCES_IO_BASE_XML} ${SOURCES_IO_BOOST_XML} ${SOURCES_IO_RAPID_XML})

if(QT4_FOUND)
    GET_SOURCE_FILES(SOURCES_IO_QT_XML IO/XmlIO/Qt)
    set(SOURCES ${SOURCES} ${SOURCES_IO_QT_XML})
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
