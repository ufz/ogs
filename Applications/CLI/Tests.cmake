
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE)

AddTest(
	NAME GroundWaterFlowProcess
	PATH Elliptic/quad_20x10_GroundWaterFlow
	EXECUTABLE ogs # optional
	EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj
	WRAPPER time # optional
	# TODO test output to expected results
	# TESTER diff
	# DIFF_DATA quad_20x10_constMat0.mesh.vtu quad_20x10_left_right.gml
	DATA quad_20x10_GroundWaterFlow.prj quad_20x10_constMat0.mesh.vtu quad_20x10_left_right.gml
)
