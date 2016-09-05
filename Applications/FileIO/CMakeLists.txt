GET_SOURCE_FILES(SOURCES)
if(NOT Shapelib_FOUND)
    list(REMOVE_ITEM SOURCES SHPInterface.h SHPInterface.cpp)
endif()
if(QT4_FOUND)
    APPEND_SOURCE_FILES(SOURCES XmlIO/Qt)
endif()
APPEND_SOURCE_FILES(SOURCES Gmsh)

if(QT4_FOUND)
    APPEND_SOURCE_FILES(SOURCES FEFLOW)
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
