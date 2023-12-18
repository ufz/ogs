add_subdirectory(ApplicationsLib)
add_subdirectory(DataHolderLib)
add_subdirectory(FileIO)

if(OGS_BUILD_UTILS AND PROJECT_IS_TOP_LEVEL)
    add_subdirectory(Utils)
elseif(OGS_BUILD_GUI)
    add_subdirectory(Utils/OGSFileConverter)
endif()

if(OGS_BUILD_GUI)
    add_subdirectory(DataExplorer)
endif() # OGS_BUILD_GUI

if(OGS_BUILD_CLI)
    add_subdirectory(CLI)
endif() # OGS_BUILD_CLI

if(OGS_USE_INSITU)
    add_subdirectory(InSituLib)
endif()

add_subdirectory(Python)
