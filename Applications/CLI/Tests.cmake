
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
		TESTER numdiff
		DIFF_DATA cube_${mesh_size}_result.dat
		DATA cube_${mesh_size}.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
	)
endforeach()

foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
	AddTest(
		NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
		PATH Elliptic/cube_1x1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS cube_${mesh_size}.prj
		WRAPPER time
		TESTER numdiff
		DIFF_DATA cube_${mesh_size}_result.dat
		DATA cube_${mesh_size}.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
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
		TESTER numdiff
		DIFF_DATA square_${mesh_size}_result.dat
		DATA square_${mesh_size}.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
	)
endforeach()

foreach(mesh_size 1e5 1e6)
	AddTest(
		NAME LARGE_GroundWaterFlowProcess_square_1x1_${mesh_size}
		PATH Elliptic/square_1x1_GroundWaterFlow
		EXECUTABLE ogs
		EXECUTABLE_ARGS square_${mesh_size}.prj
		WRAPPER time
		TESTER numdiff
		DIFF_DATA square_${mesh_size}_result.dat
		DATA square_${mesh_size}.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
	)
endforeach()

AddTest(
	NAME MaterialPropertyGroundwater
	PATH MaterialProperty/Groundwater
	EXECUTABLE ogs
	EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj
	WRAPPER time
	TESTER numdiff
	DIFF_DATA quad_20x10_result_result.dat
	DATA quad_20x10_GroundWaterFlow.prj quad_20x10_GroundWaterFlow.vtu quad_20x10_GroundWaterFlow.msh quad_20x10_left_right.gml
)
