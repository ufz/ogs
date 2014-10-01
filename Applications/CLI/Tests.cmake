
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE)

ExternalData_Add_Test(
    data
    NAME reading_GroundWaterFlow_project
    COMMAND ogs DATA{${ExternalData_SOURCE_ROOT}/Elliptic/quad_20x10_GroundWaterFlow/quad_20x10_GroundWaterFlow.prj,quad_20x10_constMat0.mesh.vtu,quad_20x10_left_right.gml}
)
