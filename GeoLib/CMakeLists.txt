# Source files
get_source_files(SOURCES)

append_source_files(SOURCES IO)

append_source_files(SOURCES IO/XmlIO/Boost)

if(OGS_BUILD_GUI)
    append_source_files(SOURCES IO/XmlIO/Qt)
endif()

# Create the library
ogs_add_library(GeoLib ${SOURCES})
ogs_add_library(GeoLibPredicates predicates.c)
set_source_files_properties(
    predicates.c
    PROPERTIES
        COMPILE_OPTIONS
        "$<$<C_COMPILER_ID:Clang,AppleClang,GNU>:-w>$<$<C_COMPILER_ID:MSVC>:/W0>"
        SKIP_UNITY_BUILD_INCLUSION TRUE
)

target_link_libraries(
    GeoLib
    PUBLIC BaseLib Eigen3::Eigen MathLib $<$<BOOL:${OGS_BUILD_GUI}>:rapidxml>
           $<$<BOOL:${OGS_BUILD_GUI}>:Qt5::Xml>
           $<$<BOOL:${OGS_BUILD_GUI}>:Qt5::XmlPatterns>
    PRIVATE range-v3 GeoLibPredicates
)

if(OGS_USE_NETCDF)
    target_link_libraries(
        GeoLib
        PRIVATE netcdf-cxx4
    )
    target_compile_definitions(GeoLib PRIVATE OGS_USE_NETCDF)
endif()

foreach(xsd OpenGeoSysGLI.xsd OpenGeoSysSTN.xsd)
    file(COPY_FILE IO/XmlIO/${xsd} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${xsd}
         ONLY_IF_DIFFERENT
    )
    install(FILES IO/XmlIO/${xsd} DESTINATION bin)
endforeach()

target_precompile_headers(GeoLib PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]] <Eigen/Core>)
