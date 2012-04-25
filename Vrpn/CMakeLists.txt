# 24/08/2010 LB Initial implementation
# CMake file for OgsVrpn

### Source files ###
SET( SOURCES
	VrpnClient.cpp
	SpaceNavigatorClient.cpp
	QSpaceNavigatorClient.cpp
	VrpnArtTrackingClient.cpp
	QVrpnArtTrackingClient.cpp
	TrackingSettingsWidget.cpp
)

### Moc header files ###
SET( MOC_HEADERS
	VrpnClient.h
	QSpaceNavigatorClient.h
	QVrpnArtTrackingClient.h
	TrackingSettingsWidget.h
)

### UI files ###
SET( UIS
	TrackingSettingsWidgetBase.ui
)

### Header files ###
SET( HEADERS
	SpaceNavigatorClient.h
	VrpnArtTrackingClient.h
)

### Qt precompiler ###
QT4_WRAP_UI( UI_HEADERS ${UIS} )
Qt4_WRAP_CPP( MOC_SOURCES ${MOC_HEADERS} )

### Include directories ###
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/Qt/Vrpn
	${CMAKE_BINARY_DIR}/Qt/Vrpn
	${VRPN_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/Qt/VtkVis
	${CMAKE_BINARY_DIR}/Qt/VtkVis
)

# Put moc files in a project folder
SOURCE_GROUP("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
SOURCE_GROUP("Moc Files" REGULAR_EXPRESSION "moc_.*")

### Create the library ###
ADD_LIBRARY( OgsVrpn STATIC
	${SOURCES}
	${HEADERS}
	${MOC_HEADERS}
	${MOC_SOURCES}
	${UIS}
)

### Link other libraries ###
#TARGET_LINK_LIBRARIES( OgsVrpn
#	${QT_LIBRARIES}
#)
