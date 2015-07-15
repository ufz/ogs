# This file initializes the required submodules
set(OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT "" CACHE STRING
	"User given submodules which should be checked out by CMake.")
if(NOT OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT)
	set(OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT "")
endif()
set(REQUIRED_SUBMODULES
	ThirdParty/quickcheck
	ThirdParty/vtkdiff
	ThirdParty/tclap
	${OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT}
)

if(OGS_BUILD_GUI)
	set(REQUIRED_SUBMODULES ${REQUIRED_SUBMODULES} ThirdParty/tetgen)
endif()

foreach(SUBMODULE ${REQUIRED_SUBMODULES})
	if(WIN32)
		set(SUBMODULE_STATE 1)
	else()
		# Check if submodule is already initialized
		execute_process(
			COMMAND ${BASH_TOOL_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/SubmoduleCheck.sh ${SUBMODULE}
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
			RESULT_VARIABLE SUBMODULE_STATE
		)
	endif()

	if(SUBMODULE_STATE EQUAL 1)
		message(STATUS "Initializing submodule ${SUBMODULE}")
		execute_process(
			COMMAND ${GIT_TOOL_PATH} submodule update --init ${SUBMODULE}
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		)
	elseif(SUBMODULE_STATE EQUAL 2)
		message(STATUS "Updating submodule ${SUBMODULE}")
		execute_process(
			COMMAND git submodule update ${SUBMODULE}
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		)
	endif()
endforeach()
