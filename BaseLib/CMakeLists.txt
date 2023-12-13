get_source_files(SOURCES)
append_source_files(SOURCES IO)
append_source_files(SOURCES IO/XmlIO)
if(OGS_BUILD_GUI)
    append_source_files(SOURCES IO/XmlIO/Qt)
endif()

# Create the library
ogs_add_library(BaseLib GENERATE_EXPORT_HEADER ${SOURCES})

target_link_libraries(
    BaseLib
    PUBLIC Boost::boost
           spdlog
           tclap
           $<$<BOOL:${MSVC}>:WinMM> # needed for timeGetTime
           $<$<BOOL:${OGS_BUILD_GUI}>:Qt5::Xml>
           $<$<BOOL:${OGS_BUILD_GUI}>:Qt5::XmlPatterns>
           $<$<BOOL:${OGS_USE_PETSC}>:PkgConfig::PETSC>
    PRIVATE xmlpatch
)

target_compile_definitions(
    BaseLib PUBLIC $<$<BOOL:${OGS_FATAL_ABORT}>:OGS_FATAL_ABORT>
)

target_precompile_headers(BaseLib PRIVATE Error.h ConfigTree.h Logging.h)
