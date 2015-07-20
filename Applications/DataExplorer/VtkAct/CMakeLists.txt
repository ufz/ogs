# Source files
set(SOURCES
	VtkPickCallback.cpp
	VtkCustomInteractorStyle.cpp
)

# Headers
set(HEADERS
	VtkPickCallback.h
	VtkCustomInteractorStyle.h
)

# Put moc files in a project folder
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")


include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/../../BaseLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../GeoLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../FileIO
	${CMAKE_CURRENT_SOURCE_DIR}/../../FileIO/GmshIO
	${CMAKE_CURRENT_SOURCE_DIR}/../../MeshLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../MathLib
	${CMAKE_CURRENT_SOURCE_DIR}/../VtkVis
)

# Create the library
add_library( VtkAct STATIC
	${SOURCES}
	${HEADERS}
)

include(AddCatalystDependency)
ADD_CATALYST_DEPENDENCY(VtkAct)

target_link_libraries( VtkAct PUBLIC Qt4::QtCore )
target_link_libraries( VtkAct INTERFACE ${CATALYST_LIBRARIES} )

set_property(TARGET VtkAct PROPERTY FOLDER "DataExplorer")
