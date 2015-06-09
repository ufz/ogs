#
# AddTest
# -------
#
# Creates application test runs. Order of arguments can be arbitrary.
#
# AddTest(
#   NAME <name of the the test>
#   PATH <working directory> # relative to SourceDir/Tests/Data
#   EXECUTABLE <executable target> # optional, defaults to ogs
#   EXECUTABLE_ARGS <arguments> # files referenced in the DATA argument can be used here
#   WRAPPER <time|memcheck|callgrind> # optional, defaults to time
#   TESTER <diff|memcheck> # optional
#   DATA <list of all required data files, white-space separated, have to be in PATH>
# )
#
# Conditional arguments:
#
#   diff-tester
#     - DIFF_DATA <list of files to diff>
#       # the given file is compared to [filename]_expected.[extension]
#
#   numdiff-tester
#     - DIFF_DATA <list of files to numdiff>
#       # the given file is compared to [filename]_expected.[extension]
#
#   vtkdiff-tester
#     - DIFF_DATA <vtk file> <data array a name> <data array b name>
#       # the given data arrays in the vtk file are compared
#

function (AddTest)

	# parse arguments
	set(options NONE)
	set(oneValueArgs EXECUTABLE PATH NAME WRAPPER TESTER)
	set(multiValueArgs EXECUTABLE_ARGS DATA DIFF_DATA WRAPPER_ARGS)
	cmake_parse_arguments(AddTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
	set(AddTest_SOURCE_PATH "${ExternalData_SOURCE_ROOT}/${AddTest_PATH}")
	set(AddTest_BINARY_PATH "${ExternalData_BINARY_ROOT}/${AddTest_PATH}")

	# set defaults
	if(NOT AddTest_EXECUTABLE)
		set(AddTest_EXECUTABLE ogs)
	endif()

	if(NOT AddTest_WRAPPER)
		set(AddTest_WRAPPER time)
	endif()

	# replace arguments which reference test data files with the correct DATA{}-path
	foreach(ARG ${AddTest_EXECUTABLE_ARGS})
		string(REGEX MATCH ".*${ARG}.*" ARG_FOUND ${AddTest_DATA} )
		if(ARG_FOUND)
			set(AddTest_EXECUTABLE_ARGS_PARSED ${AddTest_EXECUTABLE_ARGS_PARSED} DATA{${AddTest_SOURCE_PATH}/${ARG}})
		else()
			set(AddTest_EXECUTABLE_ARGS_PARSED ${AddTest_EXECUTABLE_ARGS_PARSED} ${ARG}})
		endif()
	endforeach()

	string(REPLACE ";" "," AddTest_DATA "${AddTest_DATA}")
	set(AddTest_DATA "${AddTest_SOURCE_PATH}/${AddTest_DATA}")


	# --- Implement wrappers ---
	# check requirements, disable if not met
	if(AddTest_WRAPPER STREQUAL "time" AND NOT TIME_TOOL_PATH)
		return()
	endif()
	if(AddTest_WRAPPER STREQUAL "memcheck" AND NOT VALGRIND_TOOL_PATH)
		return()
	endif()
	if(AddTest_WRAPPER STREQUAL "callgrind" AND NOT VALGRIND_TOOL_PATH)
		return()
	endif()
	if(AddTest_WRAPPER STREQUAL "mpirun" AND NOT MPIRUN_TOOL_PATH)
		return()
	endif()

	if(AddTest_WRAPPER STREQUAL "time")
		set(WRAPPER_COMMAND time)
	elseif(AddTest_WRAPPER STREQUAL "memcheck" AND VALGRIND_TOOL_PATH)
		set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=memcheck --log-file=${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.log -v --leak-check=full --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xff")
		set(tester memcheck)
	elseif(AddTest_WRAPPER STREQUAL "callgrind" AND VALGRIND_TOOL_PATH)
		set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=callgrind --branch-sim=yes --cache-sim=yes --dump-instr=yes --collect-jumps=yes")
		unset(tester)
	elseif(AddTest_WRAPPER STREQUAL "mpirun")
		set(WRAPPER_COMMAND "${MPIRUN_TOOL_PATH} ${AddTest_WRAPPER_ARGS}")
	endif()

	# --- Implement testers ---
	# check requirements, disable if not met
	if(AddTest_TESTER STREQUAL "diff" AND NOT DIFF_TOOL_PATH)
		return()
	endif()
	if(AddTest_TESTER STREQUAL "numdiff" AND NOT NUMDIFF_TOOL_PATH)
		return()
	endif()
	if(AddTest_TESTER STREQUAL "vtkdiff" AND NOT TARGET vtkdiff)
		return()
	endif()
	if(AddTest_TESTER STREQUAL "memcheck" AND NOT GREP_TOOL_PATH)
		return()
	endif()

	if((AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER STREQUAL "numdiff" OR AddTest_TESTER STREQUAL "vtkdiff") AND NOT AddTest_DIFF_DATA)
		message(FATAL_ERROR "AddTest(): ${AddTest_NAME} - no DIFF_DATA given!")
	endif()

	if(AddTest_TESTER STREQUAL "diff")
		set(SELECTED_DIFF_TOOL_PATH ${DIFF_TOOL_PATH})
		set(TESTER_ARGS "-sbB")
	elseif(AddTest_TESTER STREQUAL "numdiff")
		set(SELECTED_DIFF_TOOL_PATH ${NUMDIFF_TOOL_PATH})
		set(TESTER_ARGS "--statistics --absolute-tolerance=1e-5 --relative-tolerance=1e-4")
	elseif(AddTest_TESTER STREQUAL "vtkdiff")
		set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:vtkdiff>)
		set(TESTER_ARGS "-q --abs 1e-2 --rel 1e-4")
	endif()

	if(AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER STREQUAL "numdiff")
		foreach(FILE ${AddTest_DIFF_DATA})
			get_filename_component(FILE_NAME ${FILE} NAME_WE)
			get_filename_component(FILE_EXT ${FILE} EXT)
			set(FILE_EXPECTED ${FILE_NAME}_expected${FILE_EXT})
			set(TESTER_COMMAND ${TESTER_COMMAND} "${SELECTED_DIFF_TOOL_PATH} ${TESTER_ARGS} DATA{${AddTest_SOURCE_PATH}/${FILE_EXPECTED}} ${AddTest_BINARY_PATH}/${FILE}")
			if(AddTest_DIFF_DATA_PARSED)
				set(AddTest_DIFF_DATA_PARSED "${AddTest_DIFF_DATA_PARSED},${FILE_EXPECTED}")
			else()
				set(AddTest_DIFF_DATA_PARSED "${FILE_EXPECTED}")
			endif()
		endforeach()
		string(REPLACE ";" " && " TESTER_COMMAND "${TESTER_COMMAND}")
		set(AddTest_DIFF_DATA_PARSED "${AddTest_SOURCE_PATH}/${AddTest_DIFF_DATA_PARSED}")
	elseif(AddTest_TESTER STREQUAL "vtkdiff")
		list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        if (NOT ${DiffDataLength} EQUAL 3)
			message(FATAL_ERROR "For vtkdiff tester 3 diff data arguments are required.")
		endif()
		list(GET AddTest_DIFF_DATA 0 VTK_FILE)
		list(GET AddTest_DIFF_DATA 1 NAME_A)
		list(GET AddTest_DIFF_DATA 2 NAME_B)

		set(TESTER_COMMAND ${TESTER_COMMAND} "${SELECTED_DIFF_TOOL_PATH} \
			${AddTest_BINARY_PATH}/${VTK_FILE} -a ${NAME_A} -b ${NAME_B} \
			${TESTER_ARGS}")
		string(REPLACE ";" " && " TESTER_COMMAND "${TESTER_COMMAND}")
	elseif(tester STREQUAL "memcheck")
		set(TESTER_COMMAND "! ${GREP_TOOL_PATH} definitely ${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.log")
	endif()

	## -----------
	if(TARGET ${AddTest_EXECUTABLE})
		set(AddTest_EXECUTABLE_PARSED $<TARGET_FILE:${AddTest_EXECUTABLE}>)
	else()
		set(AddTest_EXECUTABLE_PARSED ${AddTest_EXECUTABLE})
	endif()

	# Run the wrapper
	ExternalData_Add_Test(data
		NAME "${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER}"
		COMMAND ${CMAKE_COMMAND}
		-DEXECUTABLE=${AddTest_EXECUTABLE_PARSED}
		-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS_PARSED}
		-Dcase_path=${AddTest_SOURCE_PATH}
		-Dcase_name=${AddTest_NAME}
		-DWRAPPER_COMMAND=${WRAPPER_COMMAND}
		-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
		-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
		DATA{${AddTest_DATA}}
	)

	if(NOT AddTest_TESTER)
		return()
	endif()

	# Run the tester
	if(AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER STREQUAL "numdiff")
		ExternalData_Add_Test(data
			NAME "${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER}-${AddTest_TESTER}"
			COMMAND ${CMAKE_COMMAND}
			-Dcase_path=${AddTest_SOURCE_PATH}
			-Dcase_name=${AddTest_NAME}
			-DTESTER_COMMAND=${TESTER_COMMAND}
			-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
			-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
			DATA{${AddTest_DIFF_DATA_PARSED}}
		)
	elseif(AddTest_TESTER STREQUAL "vtkdiff")
		add_test(
			NAME "${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER}-${AddTest_TESTER}"
			COMMAND ${CMAKE_COMMAND}
			-Dcase_path=${AddTest_SOURCE_PATH}
			-Dcase_name=${AddTest_NAME}
			-DTESTER_COMMAND=${TESTER_COMMAND}
			-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
			-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
		)
	else()
		add_test(
			NAME "${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER}-${AddTest_TESTER}"
			COMMAND ${CMAKE_COMMAND}
			-Dcase_path=${AddTest_SOURCE_PATH}
			-Dcase_name=${AddTest_NAME}
			-DTESTER_COMMAND=${TESTER_COMMAND}
			-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
			-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
		)
	endif()
	set_tests_properties(${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER}-${AddTest_TESTER}
		PROPERTIES DEPENDS ${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER})

endfunction()
