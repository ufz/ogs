# Source files
GET_SOURCE_FILES(SOURCES_FEMLIB)
SET ( SOURCES ${SOURCES_FEMLIB})

# Create the library
ADD_LIBRARY(FemLib STATIC ${SOURCES})

include_directories(
	.
	../BaseLib
	../GeoLib
	../MathLib
)


target_link_libraries (FemLib
	BaseLib
	GeoLib
	MathLib
)

