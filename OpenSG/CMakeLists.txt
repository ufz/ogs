# 24/08/2010 LB Initial implementation
# CMake file for OgsOpenSG

### Source files ###
SET( SOURCES
	vtkOsgActor.cpp
	OsgWidget.cpp
)

### Moc header files ###
SET( MOC_HEADERS
	OsgWidget.h
)

### UI files ###
#SET( UIS
#	widget.ui
#)

### Header files ###
SET( HEADERS
	vtkOsgActor.h
)

### Qt precompiler ###
#QT4_WRAP_UI( UI_HEADERS ${UIS} )
Qt4_WRAP_CPP( MOC_SOURCES ${MOC_HEADERS} )

### Include directories ###
INCLUDE_DIRECTORIES(
	.
	${CMAKE_BINARY_DIR}/.
)

# Put moc files in a project folder
#SOURCE_GROUP("UI Files" REGULAR_EXPRESSION "\w*\.ui")
#SOURCE_GROUP("Moc Files" REGULAR_EXPRESSION "moc_.*")

### Create the library ###
ADD_LIBRARY( OgsOpenSG STATIC
	${SOURCES}
	${HEADERS}
	${MOC_HEADERS}
	${MOC_SOURCES}
	#${UIS}
)

### Link other libraries ###
TARGET_LINK_LIBRARIES( OgsOpenSG
	${QT_LIBRARIES}
)

USE_OPENSG(OgsOpenSG)