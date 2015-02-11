add_subdirectory(ApplicationsLib)

if(OGS_BUILD_UTILS AND NOT IS_SUBPROJECT)
	add_subdirectory(Utils)
endif() # OGS_BUILD_UTILS AND NOT IS_SUBPROJECT

if(OGS_BUILD_GUI)
	add_definitions(-DOGS_BUILD_GUI)
	add_subdirectory(DataExplorer)
endif() # OGS_BUILD_GUI


if(OGS_BUILD_CLI)
    add_subdirectory(CLI)
endif() # OGS_BUILD_CLI
