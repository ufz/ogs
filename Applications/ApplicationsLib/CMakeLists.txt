# Source files
SET ( LIB_SOURCES ProjectData.cpp )

# Library
ADD_LIBRARY (ApplicationsLib STATIC ${LIB_SOURCES})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../DataExplorer/DataView)

TARGET_LINK_LIBRARIES( ApplicationsLib
	GeoLib
	ProcessLib
)
