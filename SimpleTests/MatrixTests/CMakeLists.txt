## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD On )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
        SET (HAVE_PTHREADS TRUE)
        MESSAGE (STATUS "pthread library found." )
ENDIF (CMAKE_USE_PTHREADS_INIT )

FIND_PACKAGE( Metis )
IF (METIS_FOUND)
	MESSAGE (STATUS "metis found.")
ELSE(METIS_FOUND)
	MESSAGE (STATUS "metis not found.")
ENDIF(METIS_FOUND)

INCLUDE_DIRECTORIES(
        .
	../../Base/
	../../MathLib/
)

# Create the executable
ADD_EXECUTABLE( MatMult
        MatMult.cpp
        ${SOURCES}
        ${HEADERS}
)

# Create the executable
ADD_EXECUTABLE( MatTestRemoveRowsCols
        MatTestRemoveRowsCols.cpp
        ${SOURCES}
        ${HEADERS}
)

ADD_EXECUTABLE( MatVecMultPerm
        MatVecMultPerm.cpp
        ${SOURCES}
        ${HEADERS}
)


IF (WIN32)
        TARGET_LINK_LIBRARIES(MatMult Winmm.lib)
ENDIF (WIN32)

TARGET_LINK_LIBRARIES ( MatMult
	Base
	MathLib
)

IF (WIN32)
        TARGET_LINK_LIBRARIES(MatTestRemoveRowsCols Winmm.lib)
ENDIF (WIN32)

TARGET_LINK_LIBRARIES ( MatTestRemoveRowsCols
	Base
	MathLib
)

TARGET_LINK_LIBRARIES ( MatVecMultPerm
	Base
	MathLib
	${METIS_LIBRARIES}
)
