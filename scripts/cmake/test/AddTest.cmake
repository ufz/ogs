FUNCTION (AddOgsBenchmark name project_name)

	ExternalData_Add_Test(data
		NAME "${name}-${project_name}"
		COMMAND ogs DATA{${ExternalData_SOURCE_ROOT}/${project_name}} ${ARGN}
	)

ENDFUNCTION()

FUNCTION (AddTest executable case_path case_name wrapper)

	SET(tester ${ARGV4})

	# Implement wrappers
	IF(wrapper STREQUAL "TIME")
		SET(WRAPPER_COMMAND time)
	ELSEIF(wrapper STREQUAL "MEMCHECK" AND VALGRIND_TOOL_PATH)
		SET(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=memcheck --log-file=${case_path}/${case_name}_memcheck.log -v --leak-check=full --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xff")
		SET(tester MEMCHECK)
	ELSEIF(wrapper STREQUAL "CALLGRIND" AND VALGRIND_TOOL_PATH)
		SET(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=callgrind --branch-sim=yes --cache-sim=yes --dump-instr=yes --collect-jumps=yes")
		UNSET(tester)
	ENDIF()

	# Implement testers
	IF(tester STREQUAL "DIFF")
		SET(TESTER_COMMAND "${DIFF_TOOL_PATH} -sbB ${case_path}/${case_name}_expected_result.vtu ${case_path}/${case_name}_with_results.vtu")
	ELSEIF(tester STREQUAL "MEMCHECK")
		SET(TESTER_COMMAND "! ${GREP_TOOL_PATH} definitely ${case_path}/${case_name}_memcheck.log")
	ENDIF()

	## -----------

	ExternalData_Add_Test(data
		NAME "${executable}-${case_path}-${wrapper}"
		COMMAND ${CMAKE_COMMAND}
		-Dexecutable=$<TARGET_FILE:${executable}>
		-Dcase_path=${case_path}
		-Dcase_name=${case_name}
		-Dwrapper=${wrapper}
		-DWRAPPER_COMMAND=${WRAPPER_COMMAND}
		-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
		-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
		DATA{${case_path}/${case_name}.cnd, ${case_path}/${case_name}.gml, ${case_path}/${case_name}.vtu}
	)

	IF(NOT tester)
		RETURN()
	ENDIF()

	ADD_TEST(NAME "${executable}-${case_path}-${wrapper}-${tester}"
		COMMAND ${CMAKE_COMMAND}
		-Dcase_path=${case_path}
		-Dcase_name=${case_name}
		-Dtester=${tester}
		-DTESTER_COMMAND=${TESTER_COMMAND}
		-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
		-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
	)

ENDFUNCTION()
