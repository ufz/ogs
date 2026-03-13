append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    .
    IO
    IO/XmlIO
)
if(OGS_BUILD_GUI)
    append_source_and_header_files(
        _sources _headers _private DIRECTORIES IO/XmlIO/Qt
    )
endif()

# Create the library
ogs_add_library(
    BaseLib
    ${_sources}
    GENERATE_EXPORT_HEADER
    PUBLIC_HEADERS
    ${_headers}
    PRIVATE_HEADERS
    ${_private}
)

if(MSVC)
    find_library(WinMM_LIBRARY WinMM REQUIRED)
    find_library(DbgHelp_LIBRARY DbgHelp REQUIRED)
endif()

target_link_libraries(
    BaseLib
    PUBLIC Boost::algorithm
           Boost::property_tree
           range-v3
           spdlog::spdlog
           tclap
           $<$<BOOL:${MSVC}>:${WinMM_LIBRARY}> # needed for timeGetTime
           $<$<BOOL:${MSVC}>:${DbgHelp_LIBRARY}> # needed for UnDecorateSymbolName
           $<$<BOOL:${OGS_BUILD_GUI}>:Qt5::Xml>
           $<$<BOOL:${OGS_BUILD_GUI}>:Qt5::XmlPatterns>
           $<$<BOOL:${OGS_USE_PETSC}>:PkgConfig::PETSC>
    PRIVATE xmlpatch Boost::interprocess
)

target_compile_definitions(
    BaseLib PUBLIC $<$<BOOL:${OGS_FATAL_ABORT}>:OGS_FATAL_ABORT>
)

target_precompile_headers(BaseLib PRIVATE Error.h ConfigTree.h Logging.h)
