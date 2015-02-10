## pthread ##
if (WIN32)
	set ( ADDITIONAL_LIBS Winmm.lib )
endif()

if ( CMAKE_USE_PTHREADS_INIT )
        set (HAVE_PTHREADS TRUE)
endif (CMAKE_USE_PTHREADS_INIT )

include_directories(
        .
	${CMAKE_SOURCE_DIR}/BaseLib/
	${CMAKE_BINARY_DIR}/BaseLib/
	${CMAKE_SOURCE_DIR}/MathLib/
)

# Create the executable
add_executable( MatMult
        MatMult.cpp
        ${SOURCES}
        ${HEADERS}
)
set_target_properties(MatMult PROPERTIES FOLDER SimpleTests)
target_link_libraries(MatMult logog)

if ( HAVE_PTHREADS )
add_executable( MatVecMultPthreads
        MatVecMultPthreads.cpp
        ${SOURCES}
        ${HEADERS}
)
set_target_properties(MatVecMultPthreads PROPERTIES FOLDER SimpleTests)
target_link_libraries(MatVecMultPthreads
	pthread
	BaseLib
	MathLib
	logog)
endif ( HAVE_PTHREADS )

set_target_properties(MatMult PROPERTIES FOLDER SimpleTests)
target_link_libraries(MatMult logog)

# Create the executable
add_executable( MatTestRemoveRowsCols
        MatTestRemoveRowsCols.cpp
        ${SOURCES}
        ${HEADERS}
)
set_target_properties(MatTestRemoveRowsCols PROPERTIES FOLDER SimpleTests)

if (METIS_FOUND)
	add_executable( MatVecMultNDPerm
        	MatVecMultNDPerm.cpp
	        ${SOURCES}
        	${HEADERS}
	)
	set_target_properties(MatVecMultNDPerm PROPERTIES FOLDER SimpleTests)

	target_link_libraries ( MatVecMultNDPerm
		BaseLib
		MathLib
		logog
		${METIS_LIBRARIES}
		${ADDITIONAL_LIBS}
	)

	if(OPENMP_FOUND)
		add_executable( MatVecMultNDPermOpenMP
			MatVecMultNDPermOpenMP.cpp
			${SOURCES}
			${HEADERS}
		)
		set_target_properties(MatVecMultNDPermOpenMP PROPERTIES FOLDER SimpleTests)

		target_link_libraries ( MatVecMultNDPermOpenMP
			BaseLib
			MathLib
			logog
			${METIS_LIBRARIES}
			${ADDITIONAL_LIBS}
		)
	endif()
endif(METIS_FOUND)


if (WIN32)
        target_link_libraries(MatMult Winmm.lib)
endif (WIN32)

target_link_libraries ( MatMult
	BaseLib
	MathLib
)

if (WIN32)
        target_link_libraries(MatTestRemoveRowsCols Winmm.lib)
endif (WIN32)

target_link_libraries ( MatTestRemoveRowsCols
	BaseLib
	MathLib
)

add_executable( DenseGaussEliminationChecker
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

