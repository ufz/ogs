#
# AddTest
# -------
#
# Creates application test runs. Order of arguments can be arbitrary.
#
# AddTest(
#   NAME <name of the the test>
#   PATH <working directory> # use ${ExternalData_SOURCE_ROOT}
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

function (AddTest)

	# parse arguments
	set(options NONE)
	set(oneValueArgs EXECUTABLE PATH NAME WRAPPER TESTER)
	set(multiValueArgs EXECUTABLE_ARGS DATA DIFF_DATA)
	cmake_parse_arguments(AddTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

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
			set(AddTest_EXECUTABLE_ARGS_PARSED ${AddTest_EXECUTABLE_ARGS_PARSED} DATA{${AddTest_PATH}/${ARG}})
		else()
			set(AddTest_EXECUTABLE_ARGS_PARSED ${AddTest_EXECUTABLE_ARGS_PARSED} ${ARG}})
		endif()
	endforeach()

	string(REPLACE ";" "," AddTest_DATA "${AddTest_DATA}")
	set(AddTest_DATA "${AddTest_PATH}/${AddTest_DATA}")


	# --- Implement wrappers ---
	# check requirements
	if(AddTest_WRAPPER STREQUAL "time" AND NOT TIME_TOOL_PATH)
		message(FATAL_ERROR "time-command is required for time wrapper but was not found!")
	endif()
	if(AddTest_WRAPPER STREQUAL "memcheck" AND NOT VALGRIND_TOOL_PATH)
		message(FATAL_ERROR "Valgrind is required for memcheck wrapper but was not found!")
	endif()
	if(AddTest_WRAPPER STREQUAL "callgrind" AND NOT VALGRIND_TOOL_PATH)
		message(FATAL_ERROR "Valgrind is required for callgrind wrapper but was not found!")
	endif()

	if(AddTest_WRAPPER STREQUAL "time")
		set(WRAPPER_COMMAND time)
	elseif(AddTest_WRAPPER STREQUAL "memcheck" AND VALGRIND_TOOL_PATH)
		set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=memcheck --log-file=${AddTest_PATH}/${AddTest_NAME}_memcheck.log -v --leak-check=full --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xff")
		set(tester memcheck)
	elseif(AddTest_WRAPPER STREQUAL "callgrind" AND VALGRIND_TOOL_PATH)
		set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=callgrind --branch-sim=yes --cache-sim=yes --dump-instr=yes --collect-jumps=yes")
		unset(tester)
	endif()

	# --- Implement testers ---
	# check requirements
	if(AddTest_TESTER STREQUAL "diff" AND NOT DIFF_TOOL_PATH)
		message(FATAL_ERROR "diff-command is required for diff tester but was not found!")
	endif()
	if(AddTest_TESTER STREQUAL "diff" AND NOT AddTest_DIFF_DATA)
		message(FATAL_ERROR "AddTest(): ${AddTest_NAME} - no DIFF_DATA given!")
	endif()
	if(AddTest_TESTER STREQUAL "memcheck" AND NOT GREP_TOOL_PATH)
		message(FATAL_ERROR "grep-command is required for memcheck tester but was not found!")
	endif()

	if(AddTest_TESTER STREQUAL "diff")
		foreach(FILE ${AddTest_DIFF_DATA})
			get_filename_component(FILE_NAME ${FILE} NAME_WE)
			get_filename_component(FILE_EXT ${FILE} EXT)
			set(FILE_EXPECTED ${FILE_NAME}_expected${FILE_EXT})
			set(TESTER_COMMAND ${TESTER_COMMAND} "${DIFF_TOOL_PATH} -sbB ${AddTest_PATH}/${FILE_EXPECTED} ${AddTest_PATH}/${FILE}")
			if(AddTest_DIFF_DATA_PARSED)
				set(AddTest_DIFF_DATA_PARSED "${AddTest_DIFF_DATA_PARSED},${FILE_EXPECTED}")
			else()
				set(AddTest_DIFF_DATA_PARSED "${FILE_EXPECTED}")
			endif()
		endforeach()
		string(REPLACE ";" " && " TESTER_COMMAND "${TESTER_COMMAND}")
		set(AddTest_DIFF_DATA_PARSED "${AddTest_PATH}/${AddTest_DIFF_DATA_PARSED}")
		message("foo: ${AddTest_DIFF_DATA_PARSED}")
	elseif(tester STREQUAL "memcheck")
		set(TESTER_COMMAND "! ${GREP_TOOL_PATH} definitely ${AddTest_PATH}/${AddTest_NAME}_memcheck.log")
	endif()

	## -----------

	# Run the wrapper
	ExternalData_Add_Test(data
		NAME "${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER}"
		COMMAND ${CMAKE_COMMAND}
		-DEXECUTABLE=$<TARGET_FILE:${AddTest_EXECUTABLE}>
		-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS_PARSED}
		-Dcase_path=${AddTest_PATH}
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
	if(AddTest_TESTER STREQUAL "diff")
		ExternalData_Add_Test(data
			NAME "${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER}-${AddTest_TESTER}"
			COMMAND ${CMAKE_COMMAND}
			-Dcase_path=${AddTest_PATH}
			-Dcase_name=${AddTest_NAME}
			-DTESTER_COMMAND=${TESTER_COMMAND}
			-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
			-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
			DATA{${AddTest_DIFF_DATA_PARSED}}
		)
	else()
		add_test(
			NAME "${AddTest_EXECUTABLE}-${AddTest_NAME}-${AddTest_WRAPPER}-${AddTest_TESTER}"
			COMMAND ${CMAKE_COMMAND}
			-Dcase_path=${AddTest_PATH}
			-Dcase_name=${AddTest_NAME}
			-DTESTER_COMMAND=${TESTER_COMMAND}
			-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
			-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
		)
	endif()

endfunction()
