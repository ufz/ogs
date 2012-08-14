# 24/08/2010 LB Initial implementation
# CMake file for OgsOpenSG

### Source files ###
SET( SOURCES
	vtkOsgConverter.cpp
)

### Header files ###
SET( HEADERS
	vtkOsgConverter.h
)

### Include directories ###
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/Qt/OpenSG
)

### Create the library ###
ADD_LIBRARY( OgsOpenSG STATIC
	${SOURCES}
	${HEADERS}
)

USE_OPENSG(OgsOpenSG)