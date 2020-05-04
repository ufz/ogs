get_source_files(SOURCES)
append_source_files(SOURCES IO)
append_source_files(SOURCES IO/XmlIO)
if(OGS_BUILD_GUI)
    append_source_files(SOURCES IO/XmlIO/Qt)
endif()

# Create the library
ogs_add_library(BaseLib ${SOURCES})

target_link_libraries(BaseLib PUBLIC Boost::boost spdlog::spdlog
                                     std::filesystem)

if(MSVC)
    target_link_libraries(BaseLib PUBLIC WinMM) # needed for timeGetTime
endif()

if(OGS_BUILD_GUI)
    target_link_libraries(BaseLib PUBLIC Qt5::Xml Qt5::XmlPatterns)
endif()
