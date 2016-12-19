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

if(OGS_BUILD_SWMM)
    GET_SOURCE_FILES(SOURCES_SWMM SWMM)
    SET(SOURCES ${SOURCES} ${SOURCES_SWMM})
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Create the library
add_library(ApplicationsFileIO ${SOURCES})
target_link_libraries(ApplicationsFileIO DataHolderLib)

if(QT4_FOUND)
    target_link_libraries(ApplicationsFileIO Qt5::Xml Qt5::XmlPatterns)
    if(WIN32 AND CMAKE_CROSSCOMPILING AND OPENSSL_FOUND)
        target_link_libraries(ApplicationsFileIO Qt5::Network ${OPENSSL_LIBRARIES} ws2_32)
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

include(${PROJECT_SOURCE_DIR}/scripts/cmake/packaging/InstallXmlSchemaFiles.cmake)
InstallXmlSchemaFiles("XmlIO/*.xsd")
