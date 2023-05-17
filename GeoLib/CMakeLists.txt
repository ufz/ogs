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
    PRIVATE range-v3 tet
)

if(OGS_USE_NETCDF)
    target_link_libraries(
        GeoLib
        PUBLIC ${NETCDF_LIBRARIES_CXX} ${NETCDF_LIBRARIES_C}
    )
    if(NOT OGS_USE_CONAN)
        target_include_directories(
            GeoLib SYSTEM PUBLIC ${NETCDF_INCLUDES_C} ${NETCDF_INCLUDES_CXX}
        )
    endif()
endif()

foreach(xsd OpenGeoSysGLI.xsd OpenGeoSysSTN.xsd)
    # cmake-lint: disable=E1126
    file(COPY_FILE IO/XmlIO/${xsd} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${xsd}
         ONLY_IF_DIFFERENT
    )
    install(FILES IO/XmlIO/${xsd} DESTINATION bin)
endforeach()

target_precompile_headers(GeoLib PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]] <Eigen/Core>)
