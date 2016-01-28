# This file initializes the required submodules
set(OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT "" CACHE STRING
	"User given submodules which should be checked out by CMake.")
if(NOT OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT)
	set(OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT "")
endif()
set(REQUIRED_SUBMODULES
	ThirdParty/autocheck
	ThirdParty/cmake-modules
	ThirdParty/vtkdiff
	ThirdParty/tclap
	ThirdParty/tetgen
	${OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT}
)
if(OGS_BUILD_TESTS)
	list(APPEND REQUIRED_SUBMODULES Tests/Data)
	set(Data_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Tests/Data CACHE INTERNAL "")
	set(Data_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/Tests/Data CACHE INTERNAL "")
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
