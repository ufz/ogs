## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD On )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
        SET (HAVE_PTHREADS TRUE)
        MESSAGE (STATUS "pthread library found." )
ENDIF (CMAKE_USE_PTHREADS_INIT )

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

