# Source files
SET ( LIB_SOURCES ProjectData.cpp )

# Library
ADD_LIBRARY (ApplicationsLib STATIC ${LIB_SOURCES})

TARGET_LINK_LIBRARIES( ApplicationsLib
        GeoLib
        ProcessLib
    )
