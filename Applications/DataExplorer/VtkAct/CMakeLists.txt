# Source files
set(SOURCES
	VtkPickCallback.cpp
	VtkCustomInteractorStyle.cpp
)

# Moc headers
set(MOC_HEADERS
	VtkPickCallback.h
	VtkCustomInteractorStyle.h
)

# Put moc files in a project folder
source_group("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")

# Run Qts meta object compiler moc on header files
qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_BINARY_DIR}
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
	${MOC_HEADERS}
	${MOC_SOURCES}
)

include(AddCatalystDependency)
ADD_CATALYST_DEPENDENCY(VtkAct)

target_link_libraries( VtkAct PUBLIC Qt4::QtCore )

set_property(TARGET VtkAct PROPERTY FOLDER "DataExplorer")
