if(WIN32)
	set(ADDITIONAL_LIBS Winmm.lib)
endif()

add_executable(ConjugateGradientUnpreconditioned
	ConjugateGradientUnpreconditioned.cpp
	${SOURCES}
	${HEADERS}
)
set_target_properties(ConjugateGradientUnpreconditioned PROPERTIES FOLDER SimpleTests)

add_executable(ConjugateGradientDiagPrecond
	ConjugateGradientDiagonalPreconditioned.cpp
	${SOURCES}
	${HEADERS}
)
set_target_properties(ConjugateGradientDiagPrecond PROPERTIES FOLDER SimpleTests)

add_executable(BiCGStabDiagPrecond
	BiCGStabDiagPrecond.cpp
	${SOURCES}
	${HEADERS}
)
set_target_properties(BiCGStabDiagPrecond PROPERTIES FOLDER SimpleTests)

add_executable(GMResDiagPrecond
	GMResDiagPrecond.cpp
	${SOURCES}
	${HEADERS}
)
set_target_properties(GMResDiagPrecond PROPERTIES FOLDER SimpleTests)

target_link_libraries(ConjugateGradientUnpreconditioned
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	${ADDITIONAL_LIBS}
	MathLib
	BaseLib
)

target_link_libraries(ConjugateGradientDiagPrecond
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	${ADDITIONAL_LIBS}
	MathLib
	BaseLib
)

target_link_libraries(BiCGStabDiagPrecond
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	${ADDITIONAL_LIBS}
	MathLib
	BaseLib
)

target_link_libraries(GMResDiagPrecond
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	${ADDITIONAL_LIBS}
	MathLib
	BaseLib
)
