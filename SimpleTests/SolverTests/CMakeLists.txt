## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD On )
FIND_PACKAGE( Threads )

# Find blas
FIND_PACKAGE ( BLAS REQUIRED )
# Find lapack
FIND_PACKAGE ( LAPACK REQUIRED )

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
ADD_EXECUTABLE( ConjugateGradientUnpreconditioned
        ConjugateGradientUnpreconditioned.cpp
        ${SOURCES}
        ${HEADERS}
)

ADD_EXECUTABLE( ConjugateGradientDiagPrecond
        ConjugateGradientDiagonalPreconditioned.cpp
        ${SOURCES}
        ${HEADERS}
)

ADD_EXECUTABLE( BiCGStabDiagPrecond
	BiCGStabDiagPrecond.cpp
        ${SOURCES}
        ${HEADERS}
)

TARGET_LINK_LIBRARIES ( ConjugateGradientUnpreconditioned
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	Base
)

TARGET_LINK_LIBRARIES ( ConjugateGradientDiagPrecond
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	Base
)

TARGET_LINK_LIBRARIES( BiCGStabDiagPrecond
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	Base
)

