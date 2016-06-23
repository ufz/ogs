GET_SOURCE_FILES(SOURCES_FileIO)
set(SOURCES ${SOURCES_FileIO})
if(NOT Shapelib_FOUND)
    list(REMOVE_ITEM SOURCES SHPInterface.h SHPInterface.cpp)
endif()
GET_SOURCE_FILES(SOURCES_BASE_XML XmlIO)
set(SOURCES ${SOURCES} ${SOURCES_BASE_XML})
GET_SOURCE_FILES(SOURCES_BOOST_XML XmlIO/Boost)
set(SOURCES ${SOURCES} ${SOURCES_BOOST_XML})
if(QT4_FOUND)
    GET_SOURCE_FILES(SOURCES_QT_XML XmlIO/Qt)
    set(SOURCES ${SOURCES} ${SOURCES_QT_XML})
endif()
GET_SOURCE_FILES(SOURCES_RAPID_XML XmlIO/Rapid)
set(SOURCES ${SOURCES} ${SOURCES_RAPID_XML})

if(QT4_FOUND)
    GET_SOURCE_FILES(SOURCES_FEFLOW FEFLOW)
    SET(SOURCES ${SOURCES} ${SOURCES_FEFLOW})
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Create the library
add_library(ApplicationsFileIO ${SOURCES})
target_link_libraries(ApplicationsFileIO DataHolderLib)

if(QT4_FOUND)
    target_link_libraries(ApplicationsFileIO Qt4::QtXml Qt4::QtXmlPatterns)
    if(WIN32 AND CMAKE_CROSSCOMPILING AND OPENSSL_FOUND)
        target_link_libraries(ApplicationsFileIO Qt4::QtNetwork ${OPENSSL_LIBRARIES} ws2_32)
    endif()
endif()

if(Shapelib_FOUND)
    target_link_libraries(ApplicationsFileIO ${Shapelib_LIBRARIES})
endif()

ADD_VTK_DEPENDENCY(ApplicationsFileIO)

if(TARGET Eigen)
    add_dependencies(ApplicationsFileIO Eigen)
endif()

if(TARGET Boost)
    add_dependencies(ApplicationsFileIO Boost)
endif()

file(GLOB XSD_FILES . XmlIO/*.xsd)
if(APPLE AND OGS_BUILD_GUI)
    install(FILES ${XSD_FILES} DESTINATION ${CMAKE_BINARY_DIR}/_CPack_Packages/Darwin/DragNDrop/${CPACK_PACKAGE_FILE_NAME}/ALL_IN_ONE/DataExplorer.app/Contents/MacOS COMPONENT ogs_gui)
else()
    install(FILES ${XSD_FILES} DESTINATION bin COMPONENT ogs_cli)
endif()
