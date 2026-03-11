append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    .
    Gmsh
    Legacy
    GocadIO
)

if(NOT OGS_BUILD_GUI)
    list(REMOVE_ITEM _headers SHPInterface.h)
    list(REMOVE_ITEM _sources SHPInterface.cpp)
endif()

if(OGS_BUILD_GUI)
    append_source_and_header_files(
        _sources _headers _private DIRECTORIES XmlIO/Qt FEFLOW
    )
endif()

if(OGS_BUILD_SWMM)
    append_source_and_header_files(_sources _headers _private DIRECTORIES SWMM)
endif()

# Create the library
ogs_add_library(ApplicationsFileIO ${_sources} PUBLIC_HEADERS ${_headers})
target_link_libraries(
    ApplicationsFileIO
    PUBLIC BaseLib
           DataHolderLib
           GeoLib
           MathLib
           # Needed for the XmlPrjInterface, which links the DE/Base/OGSError.h.
           $<$<BOOL:${OGS_BUILD_GUI}>:QtBase>
           $<$<TARGET_EXISTS:shp>:shp>
           $<$<TARGET_EXISTS:SwmmInterface>:SwmmInterface>
           Boost::dynamic_bitset
    PRIVATE MeshLib MeshToolsLib GitInfoLib
)

target_precompile_headers(ApplicationsFileIO PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]] [["MeshLib/Mesh.h"]])

foreach(xsd OpenGeoSysCND.xsd OpenGeoSysNum.xsd OpenGeoSysProject.xsd)
    file(COPY_FILE XmlIO/${xsd} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${xsd}
         ONLY_IF_DIFFERENT
    )
    install(FILES XmlIO/${xsd} DESTINATION bin)
endforeach()
