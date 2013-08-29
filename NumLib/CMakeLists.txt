# Source files
GET_SOURCE_FILES(SOURCES_NUMLIB)
SET ( SOURCES ${SOURCES_NUMLIB})
GET_SOURCE_FILES(SOURCES_FEM Fem)
SET ( SOURCES ${SOURCES_FEM})
GET_SOURCE_FILES(SOURCES_FEM_Element Fem/FiniteElement)
SET ( SOURCES ${SOURCES} ${SOURCES_FEM_Element})
GET_SOURCE_FILES(SOURCES_FEM_Integration Fem/Integration)
SET ( SOURCES ${SOURCES} ${SOURCES_FEM_Integration})
GET_SOURCE_FILES(SOURCES_FEM_ShapeFunction Fem/ShapeFunction)
SET ( SOURCES ${SOURCES} ${SOURCES_FEM_ShapeFunction})
GET_SOURCE_FILES(SOURCES_FEM_Extrapolation Fem/Extrapolation)
SET ( SOURCES ${SOURCES} ${SOURCES_FEM_Extrapolation})
GET_SOURCE_FILES(SOURCES_FEM_CoordinatesMapping Fem/CoordinatesMapping)
SET ( SOURCES ${SOURCES} ${SOURCES_FEM_CoordinatesMapping})

# Create the library
ADD_LIBRARY(NumLib STATIC ${SOURCES})

include_directories(
	.
	${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/BaseLib
    ${CMAKE_SOURCE_DIR}/GeoLib
    ${CMAKE_SOURCE_DIR}/MathLib
    ${CMAKE_SOURCE_DIR}/MeshLib
)

IF (OGS_USE_EIGEN)
	INCLUDE_DIRECTORIES (${EIGEN_INCLUDE_DIR})
ENDIF()

target_link_libraries (NumLib
	BaseLib
	GeoLib
	MathLib
	MeshLib
)

