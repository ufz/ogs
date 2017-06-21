GET_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES Gmsh)

if(NOT Shapelib_FOUND)
    list(REMOVE_ITEM SOURCES SHPInterface.h SHPInterface.cpp)
endif()

if(Qt5XmlPatterns_FOUND)
    APPEND_SOURCE_FILES(SOURCES XmlIO/Qt)
    APPEND_SOURCE_FILES(SOURCES FEFLOW)
endif()

if(OGS_BUILD_SWMM)
    GET_SOURCE_FILES(SOURCES_SWMM SWMM)
    SET(SOURCES ${SOURCES} ${SOURCES_SWMM})
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Create the library
add_library(ApplicationsFileIO ${SOURCES})
target_link_libraries(ApplicationsFileIO
    PUBLIC BaseLib DataHolderLib GeoLib MathLib logog
    PRIVATE MeshLib
)

if(Shapelib_FOUND)
    target_link_libraries(ApplicationsFileIO PRIVATE ${Shapelib_LIBRARIES})
endif()

if(OGS_BUILD_SWMM)
    target_link_libraries(ApplicationsFileIO PRIVATE SwmmInterface)
endif()

ADD_VTK_DEPENDENCY(ApplicationsFileIO)

if(TARGET Eigen)
    add_dependencies(ApplicationsFileIO Eigen)
endif()

if(TARGET Boost)
    add_dependencies(ApplicationsFileIO Boost)
endif()

if(OGS_USE_PCH)
    cotire(ApplicationsFileIO)
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/packaging/InstallXmlSchemaFiles.cmake)
InstallXmlSchemaFiles("XmlIO/*.xsd")
