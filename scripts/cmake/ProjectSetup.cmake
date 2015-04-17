# Set build directories
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# Logging level
if(OGS_DISABLE_LOGGING)
	set(OGS_LOG_LEVEL LOGOG_LEVEL_NONE)
endif()

if(NOT DEFINED OGS_LOG_LEVEL)
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		add_definitions(-DLOGOG_LEVEL=LOGOG_LEVEL_DEBUG)
	else()
		add_definitions(-DLOGOG_LEVEL=LOGOG_LEVEL_INFO)
	endif() # CMAKE_BUILD_TYPE = Debug
else()
	add_definitions(-DLOGOG_LEVEL=${OGS_LOG_LEVEL})
endif() # NOT DEFINED OGS_LOG_LEVEL

# Enable Visual Studio project folder grouping
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

site_name(HOSTNAME)
