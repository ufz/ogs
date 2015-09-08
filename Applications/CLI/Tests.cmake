
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE)

# CUBE 1x1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3)
	AddTest(
		NAME GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
		PATH Elliptic/cube_1x1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS cube_${mesh_size}.prj
		WRAPPER time
		TESTER vtkdiff
		DIFF_DATA cube_${mesh_size}_pcs_0_ts_1.vtu Linear_1_to_minus1 Result
		DATA cube_${mesh_size}.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
	)

	AddTest(
		NAME GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
		PATH Elliptic/cube_1x1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
		WRAPPER time
		TESTER vtkdiff
		DIFF_DATA cube_${mesh_size}_neumann_pcs_0_ts_1.vtu D1_left_front_N1_right Result
		DATA cube_${mesh_size}_neumann.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
	)
endforeach()

foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
	AddTest(
		NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
		PATH Elliptic/cube_1x1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS cube_${mesh_size}.prj
		WRAPPER time
		TESTER vtkdiff
		DIFF_DATA cube_${mesh_size}_pcs_0_ts_1.vtu Linear_1_to_minus1 Result
		DATA cube_${mesh_size}.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
	)

	AddTest(
		NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
		PATH Elliptic/cube_1x1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
		WRAPPER time
		TESTER vtkdiff
		DIFF_DATA cube_${mesh_size}_neumann_pcs_0_ts_1.vtu D1_left_front_N1_right Result
		DATA cube_${mesh_size}_neumann.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
	)
endforeach()

# SQUARE 1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
	AddTest(
		NAME GroundWaterFlowProcess_square_1x1_${mesh_size}
		PATH Elliptic/square_1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS square_${mesh_size}.prj
		WRAPPER time
		TESTER vtkdiff
		DIFF_DATA square_${mesh_size}_pcs_0_ts_1.vtu Linear_1_to_minus1 Result
		DATA square_${mesh_size}.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
	)

	AddTest(
		NAME GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
		PATH Elliptic/square_1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
		WRAPPER time
		TESTER vtkdiff
		DIFF_DATA square_${mesh_size}_neumann_pcs_0_ts_1.vtu D1_left_bottom_N1_right Result
		DATA square_${mesh_size}_neumann.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
	)
endforeach()

foreach(mesh_size 1e5 1e6)
	AddTest(
		NAME LARGE_GroundWaterFlowProcess_square_1x1_${mesh_size}
		PATH Elliptic/square_1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS square_${mesh_size}.prj
		WRAPPER time
		TESTER vtkdiff
		DIFF_DATA square_${mesh_size}_pcs_0_ts_1.vtu Linear_1_to_minus1 Result
		DATA square_${mesh_size}.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
	)

	AddTest(
		NAME LARGE_GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
		PATH Elliptic/square_1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
		WRAPPER time
		TESTER vtkdiff
		DIFF_DATA square_${mesh_size}_neumann_pcs_0_ts_1.vtu D1_left_bottom_N1_right Result
		DATA square_${mesh_size}_neumann.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
	)
endforeach()
