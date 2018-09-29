include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

add_subdirectory(ApplicationsLib)
add_subdirectory(DataHolderLib)
add_subdirectory(FileIO)

if(OGS_BUILD_UTILS AND NOT IS_SUBPROJECT)
    add_subdirectory(Utils)
elseif(OGS_BUILD_GUI)
    add_subdirectory(Utils/OGSFileConverter)
endif() # OGS_BUILD_UTILS AND NOT IS_SUBPROJECT

if(OGS_BUILD_GUI)
    add_subdirectory(DataExplorer)
endif() # OGS_BUILD_GUI

if(OGS_BUILD_CLI)
    add_subdirectory(CLI)
endif() # OGS_BUILD_CLI

if(OGS_INSITU)
    add_subdirectory(InSituLib)
endif()
