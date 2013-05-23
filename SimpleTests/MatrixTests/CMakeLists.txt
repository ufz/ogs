## pthread ##
IF (WIN32)
	SET ( ADDITIONAL_LIBS Winmm.lib )
ENDIF()

IF ( CMAKE_USE_PTHREADS_INIT )
        SET (HAVE_PTHREADS TRUE)
ENDIF (CMAKE_USE_PTHREADS_INIT )

INCLUDE_DIRECTORIES(
        .
	${CMAKE_SOURCE_DIR}/BaseLib/
	${CMAKE_BINARY_DIR}/BaseLib/
	${CMAKE_SOURCE_DIR}/MathLib/
)

# Create the executable
ADD_EXECUTABLE( MatMult
        MatMult.cpp
        ${SOURCES}
        ${HEADERS}
)
SET_TARGET_PROPERTIES(MatMult PROPERTIES FOLDER SimpleTests)
TARGET_LINK_LIBRARIES(MatMult logog)

IF ( HAVE_PTHREADS )
ADD_EXECUTABLE( MatVecMultPthreads
        MatVecMultPthreads.cpp
        ${SOURCES}
        ${HEADERS}
)
SET_TARGET_PROPERTIES(MatVecMultPthreads PROPERTIES FOLDER SimpleTests)
TARGET_LINK_LIBRARIES(MatVecMultPthreads
	pthread
	BaseLib
	MathLib
	logog)
ENDIF ( HAVE_PTHREADS )

SET_TARGET_PROPERTIES(MatMult PROPERTIES FOLDER SimpleTests)
TARGET_LINK_LIBRARIES(MatMult logog)

# Create the executable
ADD_EXECUTABLE( MatTestRemoveRowsCols
        MatTestRemoveRowsCols.cpp
        ${SOURCES}
        ${HEADERS}
)
SET_TARGET_PROPERTIES(MatTestRemoveRowsCols PROPERTIES FOLDER SimpleTests)

IF (METIS_FOUND)
	ADD_EXECUTABLE( MatVecMultNDPerm
        	MatVecMultNDPerm.cpp
	        ${SOURCES}
        	${HEADERS}
	)
	SET_TARGET_PROPERTIES(MatVecMultNDPerm PROPERTIES FOLDER SimpleTests)

	TARGET_LINK_LIBRARIES ( MatVecMultNDPerm
		BaseLib
		MathLib
		logog
		${METIS_LIBRARIES}
		${ADDITIONAL_LIBS}
	)

	IF(OPENMP_FOUND)
		ADD_EXECUTABLE( MatVecMultNDPermOpenMP
			MatVecMultNDPermOpenMP.cpp
			${SOURCES}
			${HEADERS}
		)
		SET_TARGET_PROPERTIES(MatVecMultNDPermOpenMP PROPERTIES FOLDER SimpleTests)

		TARGET_LINK_LIBRARIES ( MatVecMultNDPermOpenMP
			BaseLib
			MathLib
			logog
			${METIS_LIBRARIES}
			${ADDITIONAL_LIBS}
		)
	ENDIF()
ENDIF(METIS_FOUND)


IF (WIN32)
        TARGET_LINK_LIBRARIES(MatMult Winmm.lib)
ENDIF (WIN32)

TARGET_LINK_LIBRARIES ( MatMult
	BaseLib
	MathLib
)

IF (WIN32)
        TARGET_LINK_LIBRARIES(MatTestRemoveRowsCols Winmm.lib)
ENDIF (WIN32)

TARGET_LINK_LIBRARIES ( MatTestRemoveRowsCols
	BaseLib
	MathLib
)

