#Source files grouped by a directory
GET_SOURCE_FILES(SOURCES_MATERIALLIB)
set ( SOURCES ${SOURCES_MATERIALLIB})

GET_SOURCE_FILES(SOURCES_FLUID_DENSITY Fluid/Density)
GET_SOURCE_FILES(SOURCES ${SOURCES} ${SOURCES_FLUID_DENSITY})

GET_SOURCE_FILES(SOURCES_SOLID_DENSITY Solid/Density)
GET_SOURCE_FILES(SOURCES ${SOURCES} ${SOURCES_SOLID_DENSITY})

include_directories (
	.
#	../BaseLib
)

# Create the library
add_library( MaterialLib STATIC ${SOURCES} )

set_target_properties(MaterialLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(MaterialLib 
#	BaseLib
)
