INCLUDE_DIRECTORIES(
        .
	../../BaseLib/
	../../MathLib/
)

# Create the executable
ADD_EXECUTABLE( ConjugateGradientUnpreconditioned
        ConjugateGradientUnpreconditioned.cpp
        ${SOURCES}
        ${HEADERS}
)
SET_TARGET_PROPERTIES(ConjugateGradientUnpreconditioned PROPERTIES FOLDER SimpleTests)

ADD_EXECUTABLE( ConjugateGradientDiagPrecond
        ConjugateGradientDiagonalPreconditioned.cpp
        ${SOURCES}
        ${HEADERS}
)
SET_TARGET_PROPERTIES(ConjugateGradientDiagPrecond PROPERTIES FOLDER SimpleTests)

ADD_EXECUTABLE( BiCGStabDiagPrecond
	BiCGStabDiagPrecond.cpp
        ${SOURCES}
        ${HEADERS}
)
SET_TARGET_PROPERTIES(BiCGStabDiagPrecond PROPERTIES FOLDER SimpleTests)

ADD_EXECUTABLE( GMResDiagPrecond
	GMResDiagPrecond.cpp
        ${SOURCES}
        ${HEADERS}
)
SET_TARGET_PROPERTIES(GMResDiagPrecond PROPERTIES FOLDER SimpleTests)


IF (WIN32)
        TARGET_LINK_LIBRARIES(ConjugateGradientUnpreconditioned Winmm.lib)
ENDIF (WIN32)
TARGET_LINK_LIBRARIES ( ConjugateGradientUnpreconditioned
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	BaseLib
)

IF (WIN32)
        TARGET_LINK_LIBRARIES(ConjugateGradientDiagPrecond Winmm.lib)
ENDIF (WIN32)
TARGET_LINK_LIBRARIES ( ConjugateGradientDiagPrecond
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	BaseLib
)

IF (WIN32)
        TARGET_LINK_LIBRARIES(BiCGStabDiagPrecond Winmm.lib)
ENDIF (WIN32)
TARGET_LINK_LIBRARIES( BiCGStabDiagPrecond
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	BaseLib
)

IF (WIN32)
        TARGET_LINK_LIBRARIES(GMResDiagPrecond Winmm.lib)
ENDIF (WIN32)
TARGET_LINK_LIBRARIES( GMResDiagPrecond
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	BaseLib
)

