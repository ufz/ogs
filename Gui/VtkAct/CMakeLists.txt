# Source files
set( SOURCES
	VtkPickCallback.cpp
	VtkCustomInteractorStyle.cpp
)

# Moc headers
set( MOC_HEADERS
	VtkPickCallback.h
	VtkCustomInteractorStyle.h
)

# Header files
set( HEADERS

)

# UI files
#set ( UIS
#)

# Put moc files in a project folder
source_group("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")

# Run Qts user interface compiler uic on .ui files
#qt4_wrap_ui( UI_HEADERS ${UIS} )

# Run Qts meta object compiler moc on header files
qt4_wrap_cpp( MOC_SOURCES ${MOC_HEADERS} )

include_directories(
	${CMAKE_SOURCE_DIR}/Gui/VtkAct
	${CMAKE_BINARY_DIR}/Gui/VtkAct
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FemLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/Gui/VtkVis
)

# Create the library
add_library( VtkAct STATIC
	${SOURCES}
	${HEADERS}
	${MOC_HEADERS}
	${MOC_SOURCES}
	${UIS}
)

target_link_libraries( VtkAct ${QT_LIBRARIES} )

set_property(TARGET VtkAct PROPERTY FOLDER "DataExplorer")