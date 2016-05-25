GET_SOURCE_FILES(SOURCES_FileIO)
set(SOURCES ${SOURCES_FileIO})
if(NOT Shapelib_FOUND)
    list(REMOVE_ITEM SOURCES SHPInterface.h SHPInterface.cpp)
endif()
if(NOT QT4_FOUND)
    list(REMOVE_ITEM SOURCES FEFLOWInterface.h FEFLOWInterface.cpp)
endif()

GET_SOURCE_FILES(SOURCES_LEGACY Legacy)
GET_SOURCE_FILES(SOURCES_GMSHIO GmshIO)
GET_SOURCE_FILES(SOURCES_RAPID_XML RapidXmlIO)
GET_SOURCE_FILES(SOURCES_VTK VtkIO)
set(SOURCES ${SOURCES} ${SOURCES_LEGACY} ${SOURCES_GMSHIO} ${SOURCES_RAPID_XML} ${SOURCES_VTK})

GET_SOURCE_FILES(SOURCES_BASE_XML XmlIO)
set(SOURCES ${SOURCES} ${SOURCES_BASE_XML})
GET_SOURCE_FILES(SOURCES_BOOST_XML XmlIO/Boost)
set(SOURCES ${SOURCES} ${SOURCES_BOOST_XML})
if(QT4_FOUND)
    GET_SOURCE_FILES(SOURCES_QT_XML XmlIO/Qt)
    set(SOURCES ${SOURCES} ${SOURCES_QT_XML})
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Create the library
add_library(FileIO ${SOURCES})
target_link_libraries(FileIO DataHolderLib)

if(QT4_FOUND)
    target_link_libraries(FileIO Qt4::QtXml Qt4::QtXmlPatterns)
    if(WIN32 AND CMAKE_CROSSCOMPILING AND OPENSSL_FOUND)
        target_link_libraries(FileIO Qt4::QtNetwork ${OPENSSL_LIBRARIES} ws2_32)
    endif()
endif()

if(Shapelib_FOUND)
    target_link_libraries(FileIO ${Shapelib_LIBRARIES})
endif()

ADD_VTK_DEPENDENCY(FileIO)

if(TARGET Eigen)
    add_dependencies(FileIO Eigen)
endif()

if(TARGET Boost)
    add_dependencies(FileIO Boost)
endif()

file(GLOB XSD_FILES . *.xsd)
if(APPLE AND OGS_BUILD_GUI)
    install(FILES ${XSD_FILES} DESTINATION ${CMAKE_BINARY_DIR}/_CPack_Packages/Darwin/DragNDrop/${CPACK_PACKAGE_FILE_NAME}/ALL_IN_ONE/DataExplorer.app/Contents/MacOS COMPONENT ogs_gui)
else()
    install(FILES ${XSD_FILES} DESTINATION bin COMPONENT ogs_cli)
endif()
