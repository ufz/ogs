include_directories(
        .
	../../BaseLib/
	../../MathLib/
)

# Create the executable
add_executable( ConjugateGradientUnpreconditioned
        ConjugateGradientUnpreconditioned.cpp
        ${SOURCES}
        ${HEADERS}
)
set_target_properties(ConjugateGradientUnpreconditioned PROPERTIES FOLDER SimpleTests)

add_executable( ConjugateGradientDiagPrecond
        ConjugateGradientDiagonalPreconditioned.cpp
        ${SOURCES}
        ${HEADERS}
)
set_target_properties(ConjugateGradientDiagPrecond PROPERTIES FOLDER SimpleTests)

add_executable( BiCGStabDiagPrecond
	BiCGStabDiagPrecond.cpp
        ${SOURCES}
        ${HEADERS}
)
set_target_properties(BiCGStabDiagPrecond PROPERTIES FOLDER SimpleTests)

add_executable( GMResDiagPrecond
	GMResDiagPrecond.cpp
        ${SOURCES}
        ${HEADERS}
)
set_target_properties(GMResDiagPrecond PROPERTIES FOLDER SimpleTests)


if (WIN32)
        target_link_libraries(ConjugateGradientUnpreconditioned Winmm.lib)
endif (WIN32)
target_link_libraries ( ConjugateGradientUnpreconditioned
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	BaseLib
)

if (WIN32)
        target_link_libraries(ConjugateGradientDiagPrecond Winmm.lib)
endif (WIN32)
target_link_libraries ( ConjugateGradientDiagPrecond
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	BaseLib
)

if (WIN32)
        target_link_libraries(BiCGStabDiagPrecond Winmm.lib)
endif (WIN32)
target_link_libraries( BiCGStabDiagPrecond
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	BaseLib
)

if (WIN32)
        target_link_libraries(GMResDiagPrecond Winmm.lib)
endif (WIN32)
target_link_libraries( GMResDiagPrecond
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
	MathLib
	BaseLib
)

