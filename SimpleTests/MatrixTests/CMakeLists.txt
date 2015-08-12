## pthread ##
if(WIN32)
	set(ADDITIONAL_LIBS Winmm.lib)
endif()

if(CMAKE_USE_PTHREADS_INIT)
	set(HAVE_PTHREADS TRUE)
endif()

# Create the executable
add_executable(MatMult
	MatMult.cpp
	${SOURCES}
	${HEADERS}
)
set_target_properties(MatMult PROPERTIES FOLDER SimpleTests)
target_link_libraries(MatMult logog)

if(HAVE_PTHREADS)
	add_executable(MatVecMultPthreads
		MatVecMultPthreads.cpp
		${SOURCES}
		${HEADERS}
	)
	set_target_properties(MatVecMultPthreads PROPERTIES FOLDER SimpleTests)
	target_link_libraries(MatVecMultPthreads
		pthread
		BaseLib
		MathLib
		logog
	)
endif()

set_target_properties(MatMult PROPERTIES FOLDER SimpleTests)
target_link_libraries(MatMult
	logog
	BaseLib
	MathLib
	${ADDITIONAL_LIBS}
)

add_executable(MatTestRemoveRowsCols
	MatTestRemoveRowsCols.cpp
	${SOURCES}
	${HEADERS}
)
set_target_properties(MatTestRemoveRowsCols PROPERTIES FOLDER SimpleTests)
target_link_libraries(MatTestRemoveRowsCols
	BaseLib
	MathLib
	${ADDITIONAL_LIBS}
)

if(METIS_FOUND)
	add_executable(MatVecMultNDPerm
		MatVecMultNDPerm.cpp
		${SOURCES}
		${HEADERS}
	)
	set_target_properties(MatVecMultNDPerm PROPERTIES FOLDER SimpleTests)

	target_link_libraries(MatVecMultNDPerm
		BaseLib
		MathLib
		logog
		${METIS_LIBRARIES}
		${ADDITIONAL_LIBS}
	)

	if(OPENMP_FOUND)
		add_executable(MatVecMultNDPermOpenMP
			MatVecMultNDPermOpenMP.cpp
			${SOURCES}
			${HEADERS}
		)
		set_target_properties(MatVecMultNDPermOpenMP PROPERTIES FOLDER SimpleTests)

		target_link_libraries(MatVecMultNDPermOpenMP
			BaseLib
			MathLib
			logog
			${METIS_LIBRARIES}
			${ADDITIONAL_LIBS}
		)
	endif()
endif()

add_executable(DenseGaussEliminationChecker
	DenseGaussEliminationChecker.cpp
	${SOURCES}
	${HEADERS}
)
set_target_properties(DenseGaussEliminationChecker PROPERTIES FOLDER SimpleTests)
target_link_libraries(DenseGaussEliminationChecker
	logog
	BaseLib
	MathLib
)

