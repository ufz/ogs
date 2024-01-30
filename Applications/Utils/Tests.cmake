AddTest(
    NAME MapGeometryToMeshSurface_Ammer
    PATH MeshGeoToolsLib/Ammer
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Ammer
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m Ammer-Homogen100m-Final-TopSurface.vtu -i Ammer-Rivers.gml -a -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Ammer/Ammer-Rivers-Mapped.gml
    TESTER diff
    DIFF_DATA Ammer-Rivers-Mapped.gml
)

AddTest(
    NAME MapGeometryToMeshSurface_Bode
    PATH MeshGeoToolsLib/Bode
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Bode
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m BodeComplex.msh -i BodeEZG_Fliessgewaesser.gml -a -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Bode/BodeEZG_Fliessgewaesser-Mapped.gml
    RUNTIME 7
    TESTER gmldiff
    DIFF_DATA BodeEZG_Fliessgewaesser-Mapped.gml 1e-10 1e-10
)

AddTest(
    NAME MapGeometryToMeshSurface_Naegelstedt
    PATH MeshGeoToolsLib/Naegelstedt
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Naegelstedt
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m SmallTest.vtu -i RiverNetwork.gml -a -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Naegelstedt/RiverNetwork-Mapped.gml
    RUNTIME 7
    TESTER diff
    DIFF_DATA RiverNetwork-Mapped.gml
)

AddTest(
    NAME createTetgenSmeshFromRasters
    PATH Utils/createTetgenSmeshFromRasters
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Utils/createTetgenSmeshFromRasters
    EXECUTABLE createTetgenSmeshFromRasters
    EXECUTABLE_ARGS -i sfc_mesh_9k.vtu -o ${Data_BINARY_DIR}/Utils/createTetgenSmeshFromRasters/CTSMFL-test -r rasterlist.txt
    RUNTIME 1
    TESTER numdiff
    DIFF_DATA CTSMFL-test.smesh 0 5e-13
)

AddTest(
    NAME createTetgenSmeshFromRasters-fail
    PATH Utils/createTetgenSmeshFromRasters
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Utils/createTetgenSmeshFromRasters
    EXECUTABLE createTetgenSmeshFromRasters
    EXECUTABLE_ARGS -i 3D.vtu -o ${Data_BINARY_DIR}/Utils/createTetgenSmeshFromRasters/CTSMFL-fail -r rasterlist.txt
    RUNTIME 1
    PROPERTIES WILL_FAIL true
)

AddTest(
    NAME postLIE
    PATH LIE/PostProcessing
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/LIE/PostProcessing
    EXECUTABLE postLIE
    EXECUTABLE_ARGS -i single_joint.pvd -o ${Data_BINARY_DIR}/LIE/PostProcessing/post_single_joint.pvd
    REQUIREMENTS OGS_BUILD_PROCESS_LIE
    TESTER vtkdiff
    DIFF_DATA
    expected_post_single_joint_ts_1_t_1.000000.vtu post_single_joint_ts_1_t_1.000000.vtu u u 1e-14 1e-14
)

AddTest(
    NAME postLIE3D
    PATH LIE/PostProcessing
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/LIE/PostProcessing
    EXECUTABLE postLIE
    EXECUTABLE_ARGS -i single_joint_3D.pvd -o ${Data_BINARY_DIR}/LIE/PostProcessing/post_single_joint_3D.pvd
    REQUIREMENTS OGS_BUILD_PROCESS_LIE
    TESTER vtkdiff
    DIFF_DATA
    post_single_joint_3D_ts_1_t_1.000000.vtu post_single_joint_3D_ts_1_t_1.000000.vtu u u 1e-14 1e-14
)

AddTest(
    NAME identifySubdomains_2D_Create
    PATH MeshGeoToolsLib/IdentifySubdomains
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/<PATH>
    EXECUTABLE identifySubdomains
    EXECUTABLE_ARGS -m 2D_mesh.vtu -o ${Data_BINARY_DIR}/<PATH>/new_ -- 2D_mesh_top_boundary.vtu 2D_mesh_bottom_boundary.vtu
    TESTER vtkdiff
    DIFF_DATA
    2D_mesh_top.vtu new_2D_mesh_top_boundary.vtu bulk_node_ids bulk_node_ids 0 0
    2D_mesh_top.vtu new_2D_mesh_top_boundary.vtu bulk_element_ids bulk_element_ids 0 0
    2D_mesh_bottom.vtu new_2D_mesh_bottom_boundary.vtu bulk_node_ids bulk_node_ids 0 0
    2D_mesh_bottom.vtu new_2D_mesh_bottom_boundary.vtu bulk_element_ids bulk_element_ids 0 0
)

AddTest(
    NAME identifySubdomains_2D_Check
    PATH MeshGeoToolsLib/IdentifySubdomains
    WORKING_DIRECTORY <SOURCE_PATH>
    EXECUTABLE identifySubdomains
    EXECUTABLE_ARGS -m 2D_mesh.vtu -o <BUILD_PATH>/check_ -- 2D_mesh_top.vtu 2D_mesh_bottom.vtu
    TESTER vtkdiff
    DIFF_DATA
    2D_mesh_top.vtu check_2D_mesh_top.vtu bulk_node_ids bulk_node_ids 0 0
    2D_mesh_top.vtu check_2D_mesh_top.vtu bulk_element_ids bulk_element_ids 0 0
    2D_mesh_bottom.vtu check_2D_mesh_bottom.vtu bulk_node_ids bulk_node_ids 0 0
    2D_mesh_bottom.vtu check_2D_mesh_bottom.vtu bulk_element_ids bulk_element_ids 0 0
)

AddTest(
    NAME identifySubdomains_riverTriangleMesh
    PATH MeshGeoToolsLib/IdentifySubdomains
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/IdentifySubdomains
    EXECUTABLE identifySubdomains
    EXECUTABLE_ARGS -m river_domain_triangle.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/IdentifySubdomains/triangle_ -- river_bc.vtu
    TESTER vtkdiff
    DIFF_DATA
    river_bc_triangle.vtu triangle_river_bc.vtu bulk_node_ids bulk_node_ids 0 0
    #river_bc_triangle.vtu triangle_river_bc.vtu bulk_element_ids bulk_element_ids 0 0   # TODO (naumov) Needs extension of vtkdiff to FieldData
    river_bc_triangle.vtu triangle_river_bc.vtu number_bulk_elements number_bulk_elements 0 0
)

AddTest(
    NAME identifySubdomains_riverPrismMesh
    PATH MeshGeoToolsLib/IdentifySubdomains
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/IdentifySubdomains
    EXECUTABLE identifySubdomains
    EXECUTABLE_ARGS -s 1e-3 -m river_domain_prism.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/IdentifySubdomains/prism_ -- river_bc.vtu
    TESTER vtkdiff
    DIFF_DATA
    river_bc_prism.vtu prism_river_bc.vtu bulk_node_ids bulk_node_ids 0 0
    #river_bc_prism.vtu prism_river_bc.vtu bulk_element_ids bulk_element_ids 0 0 # TODO (naumov) Needs extension of vtkdiff to FieldData
    river_bc_prism.vtu prism_river_bc.vtu number_bulk_elements number_bulk_elements 0 0
)

AddTest(
    NAME partmesh_2Dmesh_ogs2metis
    PATH NodePartitionedMesh/partmesh_2Dmesh_3partitions/Binary
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/partmesh_2Dmesh_3partitions/Binary
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i 2Dmesh.vtu --ogs2metis
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/partmesh_2Dmesh_3partitions/Binary
)

AddTest(
    NAME partmesh_2Dmesh_3partitions_binary
    PATH NodePartitionedMesh/partmesh_2Dmesh_3partitions/Binary
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/partmesh_2Dmesh_3partitions/Binary
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 3 -i 2Dmesh.vtu
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/partmesh_2Dmesh_3partitions/Binary --
                    2Dmesh_PLY_EAST.vtu
                    2Dmesh_PLY_WEST.vtu
                    2Dmesh_PLY_NORTH.vtu
                    2Dmesh_PLY_SOUTH.vtu
                    2Dmesh_POINT4.vtu
                    2Dmesh_POINT5.vtu
    # Mac is producing slightly different partitioning, so the results are not
    # comparable.
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_2Dmesh_ogs2metis
    TESTER diff
    DIFF_DATA 2Dmesh_partitioned_node_properties_val3.bin
              2Dmesh_partitioned_node_properties_cfg3.bin
              2Dmesh_partitioned_msh_cfg3.bin
              2Dmesh_partitioned_cell_properties_val3.bin
              2Dmesh_partitioned_cell_properties_cfg3.bin
              2Dmesh_partitioned_msh_ele_g3.bin
              2Dmesh_partitioned_msh_ele3.bin
              2Dmesh_partitioned_msh_nod3.bin
              2Dmesh_PLY_EAST_partitioned_msh_cfg3.bin
              2Dmesh_PLY_EAST_partitioned_msh_ele3.bin
              2Dmesh_PLY_EAST_partitioned_msh_ele_g3.bin
              2Dmesh_PLY_EAST_partitioned_msh_nod3.bin
              2Dmesh_PLY_NORTH_partitioned_msh_cfg3.bin
              2Dmesh_PLY_NORTH_partitioned_msh_ele3.bin
              #2Dmesh_PLY_NORTH_partitioned_msh_ele_g3.bin   empty
              2Dmesh_PLY_NORTH_partitioned_msh_nod3.bin
              2Dmesh_PLY_SOUTH_partitioned_msh_cfg3.bin
              2Dmesh_PLY_SOUTH_partitioned_msh_ele3.bin
              #2Dmesh_PLY_SOUTH_partitioned_msh_ele_g3.bin   empty
              2Dmesh_PLY_SOUTH_partitioned_msh_nod3.bin
              2Dmesh_PLY_WEST_partitioned_msh_cfg3.bin
              2Dmesh_PLY_WEST_partitioned_msh_ele3.bin
              2Dmesh_PLY_WEST_partitioned_msh_ele_g3.bin
              2Dmesh_PLY_WEST_partitioned_msh_nod3.bin
              2Dmesh_POINT4_partitioned_msh_cfg3.bin
              2Dmesh_POINT4_partitioned_msh_ele3.bin
              #2Dmesh_PLY_POINT4_partitioned_msh_ele_g3.bin   empty
              2Dmesh_POINT4_partitioned_msh_nod3.bin
              2Dmesh_POINT5_partitioned_msh_cfg3.bin
              2Dmesh_POINT5_partitioned_msh_ele3.bin
              #2Dmesh_PLY_POINT5_partitioned_msh_ele_g3.bin   empty
              2Dmesh_POINT5_partitioned_msh_nod3.bin

              2Dmesh_PLY_EAST_partitioned_cell_properties_cfg3.bin
              2Dmesh_PLY_EAST_partitioned_cell_properties_val3.bin
              2Dmesh_PLY_NORTH_partitioned_cell_properties_cfg3.bin
              2Dmesh_PLY_NORTH_partitioned_cell_properties_val3.bin
              2Dmesh_PLY_SOUTH_partitioned_cell_properties_cfg3.bin
              2Dmesh_PLY_SOUTH_partitioned_cell_properties_val3.bin
              2Dmesh_PLY_WEST_partitioned_cell_properties_cfg3.bin
              2Dmesh_PLY_WEST_partitioned_cell_properties_val3.bin
              2Dmesh_POINT4_partitioned_cell_properties_cfg3.bin
              2Dmesh_POINT4_partitioned_cell_properties_val3.bin
              2Dmesh_POINT5_partitioned_cell_properties_cfg3.bin
              2Dmesh_POINT5_partitioned_cell_properties_val3.bin

              2Dmesh_PLY_EAST_partitioned_node_properties_cfg3.bin
              2Dmesh_PLY_EAST_partitioned_node_properties_val3.bin
              2Dmesh_PLY_NORTH_partitioned_node_properties_cfg3.bin
              2Dmesh_PLY_NORTH_partitioned_node_properties_val3.bin
              2Dmesh_PLY_SOUTH_partitioned_node_properties_cfg3.bin
              2Dmesh_PLY_SOUTH_partitioned_node_properties_val3.bin
              2Dmesh_PLY_WEST_partitioned_node_properties_cfg3.bin
              2Dmesh_PLY_WEST_partitioned_node_properties_val3.bin
              2Dmesh_POINT4_partitioned_node_properties_cfg3.bin
              2Dmesh_POINT4_partitioned_node_properties_val3.bin
              2Dmesh_POINT5_partitioned_node_properties_cfg3.bin
              2Dmesh_POINT5_partitioned_node_properties_val3.bin
)

AddTest(
    NAME partmesh_mesh_for_QuadraticElements_quad8_ogs2metis
    PATH NodePartitionedMesh/QuadraticElements/Quad8
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad8
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i quad_q.vtu --ogs2metis -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad8
)

AddTest(
    NAME partmesh_mesh_for_QuadraticElements_quad8
    PATH NodePartitionedMesh/QuadraticElements/Quad8
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad8
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m  -n 2 -i quad_q.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad8
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_for_QuadraticElements_quad8_ogs2metis
    TESTER diff
    DIFF_DATA quad_q_partitioned_msh_ele2.bin
              quad_q_partitioned_msh_ele_g2.bin
              quad_q_partitioned_msh_nod2.bin
              quad_q_partitioned_msh_cfg2.bin
              quad_q_partitioned_cell_properties_cfg2.bin
              quad_q_partitioned_cell_properties_val2.bin
)

AddTest(
    NAME partmesh_mesh_for_QuadraticElements_quad9_ogs2metis
    PATH NodePartitionedMesh/QuadraticElements/Quad9
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad9
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i quad_9node.vtu --ogs2metis -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad9
)

AddTest(
    NAME partmesh_mesh_for_QuadraticElements_quad9
    PATH NodePartitionedMesh/QuadraticElements/Quad9
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad9
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 2 -i quad_9node.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad9
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_for_QuadraticElements_quad9_ogs2metis
    TESTER diff
    DIFF_DATA  quad_9node_partitioned_msh_ele2.bin
               quad_9node_partitioned_msh_ele_g2.bin
               quad_9node_partitioned_msh_nod2.bin
               quad_9node_partitioned_msh_cfg2.bin
               quad_9node_partitioned_cell_properties_cfg2.bin
               quad_9node_partitioned_cell_properties_val2.bin
)

##############Quadratic Triangle##############
AddTest(
    NAME partmesh_mesh_for_QuadraticTriangle_ogsmetis
    PATH NodePartitionedMesh/QuadraticElements/Quad_triangle
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad_triangle
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i basicQuadTri.vtu --ogs2metis -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad_triangle
)

AddTest(
    NAME partmesh_mesh_for_QuadraticTriangle
    PATH NodePartitionedMesh/QuadraticElements/Quad_triangle
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad_triangle
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 2 -i basicQuadTri.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad_triangle
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_for_QuadraticTriangle_ogsmetis
    TESTER diff
    DIFF_DATA  basicQuadTri_partitioned_cell_properties_cfg2.bin
               basicQuadTri_partitioned_msh_cfg2.bin
               basicQuadTri_partitioned_msh_ele_g2.bin
               basicQuadTri_partitioned_cell_properties_val2.bin
               basicQuadTri_partitioned_msh_ele2.bin
               basicQuadTri_partitioned_msh_nod2.bin
)
################################################

##############Quadratic Tet#####################
AddTest(
    NAME partmesh_mesh_for_QuadraticTet_ogsmetis
    PATH NodePartitionedMesh/QuadraticElements/Quad_tet
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad_tet
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i basicQuadTet.vtu --ogs2metis -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad_tet
)

AddTest(
    NAME partmesh_mesh_for_QuadraticTet
    PATH NodePartitionedMesh/QuadraticElements/Quad_tet
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad_tet
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 2 -i basicQuadTet.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad_tet
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_for_QuadraticTet_ogsmetis
    TESTER diff
    DIFF_DATA   basicQuadTet_partitioned_cell_properties_cfg2.bin
                basicQuadTet_partitioned_cell_properties_val2.bin
                basicQuadTet_partitioned_msh_cfg2.bin
                basicQuadTet_partitioned_msh_ele2.bin
                basicQuadTet_partitioned_msh_ele_g2.bin
                basicQuadTet_partitioned_msh_nod2.bin
)
################################################

##############Quadratic Hex#####################
AddTest(
    NAME partmesh_mesh_for_QuadraticHex_ogsmetis
    PATH NodePartitionedMesh/QuadraticElements/Quad_hex
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad_hex
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i basicQuadHex.vtu --ogs2metis -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad_hex
)

AddTest(
    NAME partmesh_mesh_for_QuadraticHex
    PATH NodePartitionedMesh/QuadraticElements/Quad_hex
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad_hex
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 2 -i basicQuadHex.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad_hex
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_for_QuadraticHex_ogsmetis
    TESTER diff
    DIFF_DATA   basicQuadHex_partitioned_cell_properties_cfg2.bin
                basicQuadHex_partitioned_cell_properties_val2.bin
                basicQuadHex_partitioned_msh_cfg2.bin
                basicQuadHex_partitioned_msh_ele2.bin
                basicQuadHex_partitioned_msh_ele_g2.bin
                basicQuadHex_partitioned_msh_nod2.bin
)
################################################

##############Quadratic Line#####################
AddTest(
    NAME partmesh_mesh_for_QuadraticLine_ogsmetis
    PATH NodePartitionedMesh/QuadraticElements/Quad_line
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad_line
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i basicQuadLine.vtu --ogs2metis -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad_line
)

AddTest(
    NAME partmesh_mesh_for_QuadraticLine
    PATH NodePartitionedMesh/QuadraticElements/Quad_line
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/QuadraticElements/Quad_line
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 2 -i basicQuadLine.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/QuadraticElements/Quad_line
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_for_QuadraticLine_ogsmetis
    TESTER diff
    DIFF_DATA   basicQuadLine_partitioned_cell_properties_cfg2.bin
                basicQuadLine_partitioned_cell_properties_val2.bin
                basicQuadLine_partitioned_msh_cfg2.bin
                basicQuadLine_partitioned_msh_ele2.bin
                basicQuadLine_partitioned_msh_ele_g2.bin
                basicQuadLine_partitioned_msh_nod2.bin
)

# Mesh with integration point data
AddTest(
    NAME partmesh_mesh_with_sigm_ip_pointheatsource_quad_ogs2metis
    PATH NodePartitionedMesh/WithIntegrationPointStress
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000.vtu --ogs2metis
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress
)

AddTest(
    NAME partmesh_mesh_with_sigm_ip_pointheatsource_quad
    PATH NodePartitionedMesh/WithIntegrationPointStress
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 4 -i expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_with_sigm_ip_pointheatsource_quad_ogs2metis
    TESTER diff
    DIFF_DATA
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_cell_properties_cfg4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_cell_properties_val4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_integration_point_properties_cfg4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_integration_point_properties_val4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_msh_cfg4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_msh_ele4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_msh_ele_g4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_msh_nod4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_node_properties_cfg4.bin
        expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000_partitioned_node_properties_val4.bin
)

AddTest(
    NAME partmesh_mesh_withsigma_ip_for_m1_3Dload_hex_ogs2metis
    PATH NodePartitionedMesh/WithIntegrationPointStress
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i m1_3Dload_ts_1_t_1_000000.vtu --ogs2metis
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress
)

AddTest(
    NAME partmesh_mesh_withsigma_ip_for_m1_3Dload_hex
    PATH NodePartitionedMesh/WithIntegrationPointStress
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 3 -i m1_3Dload_ts_1_t_1_000000.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_withsigma_ip_for_m1_3Dload_hex_ogs2metis
    TESTER diff
    DIFF_DATA
        m1_3Dload_ts_1_t_1_000000_partitioned_cell_properties_cfg3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_cell_properties_val3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_integration_point_properties_cfg3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_integration_point_properties_val3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_msh_cfg3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_msh_ele3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_msh_ele_g3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_msh_nod3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_node_properties_cfg3.bin
        m1_3Dload_ts_1_t_1_000000_partitioned_node_properties_val3.bin
)

AddTest(
    NAME partmesh_mesh_withsigma_ip_with_mixed_element_types_ogs2metis
    PATH NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i mesh_with_3D_different_elements_sigma_ip.vtu --ogs2metis
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements
)

AddTest(
    NAME partmesh_mesh_withsigma_ip_with_mixed_element_types
    PATH NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 2 -i mesh_with_3D_different_elements_sigma_ip.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_withsigma_ip_with_mixed_element_types_ogs2metis
    TESTER diff
    DIFF_DATA
       mesh_with_3D_different_elements_sigma_ip_partitioned_cell_properties_cfg2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_cell_properties_val2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_integration_point_properties_cfg2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_integration_point_properties_val2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_msh_cfg2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_msh_ele2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_msh_ele_g2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_msh_nod2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_node_properties_cfg2.bin
       mesh_with_3D_different_elements_sigma_ip_partitioned_node_properties_val2.bin
)

AddTest(
    NAME partmesh_mesh_with_sigm_ip_quad_tri_ogs2metis
    PATH NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i quad_tri_THM_t_864000_000000.vtu --ogs2metis
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
)

AddTest(
    NAME partmesh_mesh_with_sigm_ip_quad_tri
    PATH NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 4 -i quad_tri_THM_t_864000_000000.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_with_sigm_ip_quad_tri_ogs2metis
    TESTER diff
    DIFF_DATA
       quad_tri_THM_t_864000_000000_partitioned_cell_properties_cfg4.bin
       quad_tri_THM_t_864000_000000_partitioned_cell_properties_val4.bin
       quad_tri_THM_t_864000_000000_partitioned_integration_point_properties_cfg4.bin
       quad_tri_THM_t_864000_000000_partitioned_integration_point_properties_val4.bin
       quad_tri_THM_t_864000_000000_partitioned_msh_cfg4.bin
       quad_tri_THM_t_864000_000000_partitioned_msh_ele4.bin
       quad_tri_THM_t_864000_000000_partitioned_msh_ele_g4.bin
       quad_tri_THM_t_864000_000000_partitioned_msh_nod4.bin
       quad_tri_THM_t_864000_000000_partitioned_node_properties_cfg4.bin
       quad_tri_THM_t_864000_000000_partitioned_node_properties_val4.bin
)

AddTest(
    NAME partmesh_mesh_with_field_data_without_ip_data_ogs2metis
    PATH NodePartitionedMesh/FieldDataWithoutIPData
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/FieldDataWithoutIPData
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -i A2_tunnel_surface.vtu --ogs2metis
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/FieldDataWithoutIPData
)

AddTest(
    NAME partmesh_mesh_with_field_data_without_ip_data
    PATH NodePartitionedMesh/FieldDataWithoutIPData
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/FieldDataWithoutIPData
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 2 -i A2_tunnel_surface.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/FieldDataWithoutIPData
    REQUIREMENTS NOT APPLE
    DEPENDS partmesh-partmesh_mesh_with_field_data_without_ip_data_ogs2metis
    TESTER diff
    DIFF_DATA
        A2_tunnel_surface_partitioned_cell_properties_cfg2.bin
        A2_tunnel_surface_partitioned_msh_cfg2.bin
        A2_tunnel_surface_partitioned_node_properties_cfg2.bin
        A2_tunnel_surface_partitioned_cell_properties_val2.bin
        A2_tunnel_surface_partitioned_msh_ele2.bin
        A2_tunnel_surface_partitioned_node_properties_val2.bin
        A2_tunnel_surface_partitioned_integration_point_properties_cfg2.bin
        A2_tunnel_surface_partitioned_msh_ele_g2.bin
        A2_tunnel_surface_partitioned_integration_point_properties_val2.bin
        A2_tunnel_surface_partitioned_msh_nod2.bin
)

################################################

if(SNAKEMAKE AND TEE_TOOL_PATH AND BASH_TOOL_PATH)
    add_test(NAME snakemake_partmesh_mixed_elements
        COMMAND bash -c "export PATH=$<TARGET_FILE_DIR:partmesh>:$PATH && ${SNAKEMAKE} -j 4 \
            --config input_dir=${Data_SOURCE_DIR}/Utils/GMSH2OGS \
            -s ${PROJECT_SOURCE_DIR}/scripts/snakemake/workflows/partmesh.smk \
            ${Data_BINARY_DIR}/Utils/GMSH2OGS/{linear,quadratic}_mesh/{2,4,8,12}"
    )
    set_tests_properties(snakemake_partmesh_mixed_elements
        PROPERTIES LABELS "default")
endif()

if(SNAKEMAKE AND TEE_TOOL_PATH AND BASH_TOOL_PATH AND OGS_USE_MPI)
    add_test(NAME snakemake_reorder_mesh
        COMMAND bash -c "${SNAKEMAKE} -j 4 \
            --configfile ${PROJECT_BINARY_DIR}/buildinfo.yaml --forceall \
            -s ${PROJECT_SOURCE_DIR}/Applications/Utils/TestReorderMesh.smk"
    )
    set_tests_properties(snakemake_reorder_mesh PROPERTIES
        LABELS "default"
        RUN_SERIAL TRUE # prevent unidentified race condition
    )
endif()

# Regression test for https://gitlab.opengeosys.org/ogs/ogs/-/issues/1845 fixed in
# https://github.com/ufz/ogs/pull/2237
# checkMesh crashed when encountered Line3 element.
AddTest(
    NAME checkMesh_LIE_HM_TaskB
    PATH LIE/HydroMechanics
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/LIE/HydroMechanics
    EXECUTABLE checkMesh
    EXECUTABLE_ARGS -p -v TaskB_mesh.vtu
)

AddTest(
    NAME Mesh2Raster_small_Test
    PATH FileIO/Mesh2Raster
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileIO/Mesh2Raster
    EXECUTABLE Mesh2Raster
    EXECUTABLE_ARGS -i tri_8_ascii.vtu -o ${Data_BINARY_DIR}/FileIO/Mesh2Raster/tri_8_ascii_raster.asc -c 0.1
    TESTER diff
    DIFF_DATA tri_8_ascii_raster.asc
)

AddTest(
    NAME Mesh2Raster_large_Test
    PATH MeshGeoToolsLib/Hamburg
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Hamburg
    EXECUTABLE Mesh2Raster
    EXECUTABLE_ARGS -i 00-surface.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Hamburg/00-raster.asc -c 25
    TESTER diff
    DIFF_DATA 00-raster.asc
)

AddTest(
    NAME ExtractSurfaceLeft
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3_layers_10.vtu -o ${Data_BINARY_DIR}/MeshLib/Left.vtu -x 1 -y 0 -z 0 -a 25
    TESTER vtkdiff-mesh
    DIFF_DATA Left.vtu Left.vtu 1e-16
)

AddTest(
    NAME ExtractSurfaceRight
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3_layers_10.vtu -o ${Data_BINARY_DIR}/MeshLib/Right.vtu -x -1 -y 0 -z 0 -a 25
    TESTER vtkdiff-mesh
    DIFF_DATA Right.vtu Right.vtu 1e-16
)

AddTest(
    NAME ExtractSurfaceFront
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3_layers_10.vtu -o ${Data_BINARY_DIR}/MeshLib/Front.vtu -x 0 -y 1 -z 0 -a 25
    TESTER vtkdiff-mesh
    DIFF_DATA Front.vtu Front.vtu 1e-16
)

AddTest(
    NAME ExtractSurfaceBack
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3_layers_10.vtu -o ${Data_BINARY_DIR}/MeshLib/Back.vtu -x 0 -y -1 -z 0 -a 25
    TESTER vtkdiff-mesh
    DIFF_DATA Back.vtu Back.vtu 1e-16
)

AddTest(
    NAME ExtractSurface_QuadraticElement_Top
    PATH Utils/ExtractSurface
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Utils/ExtractSurface
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i ${Data_SOURCE_DIR}/Utils/GMSH2OGS/quadratic_mesh.vtu -o ${Data_BINARY_DIR}/Utils/ExtractSurface/quadratic_mesh_top_surface.vtu -x 0 -y 0 -z 1
    TESTER vtkdiff-mesh
    DIFF_DATA quadratic_mesh_top_surface.vtu quadratic_mesh_top_surface.vtu 1e-16
)

AddTest(
    NAME ExtractSurface_QuadraticElement_Bottom
    PATH Utils/ExtractSurface
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Utils/ExtractSurface
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i ${Data_SOURCE_DIR}/Utils/GMSH2OGS/quadratic_mesh.vtu -o ${Data_BINARY_DIR}/Utils/ExtractSurface/quadratic_mesh_bottom_surface.vtu -x 0 -y 0 -z -1
    TESTER vtkdiff-mesh
    DIFF_DATA quadratic_mesh_bottom_surface.vtu quadratic_mesh_bottom_surface.vtu 1e-16
)

AddTest(
    NAME GocadTSurface_Mesh_Test
    PATH MeshLib
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE GocadTSurfaceReader
    EXECUTABLE_ARGS -i Top-Lower-Shaly.ts -o ${Data_BINARY_DIR}/MeshLib -b
    TESTER vtkdiff-mesh
    DIFF_DATA Top-Lower-Shaly.vtu Top-Lower-Shaly.vtu 1e-16
)

if(TEST GocadTSurfaceReader-GocadTSurface_Mesh_Test-vtkdiff-mesh)
    AddTest(
        NAME GocadTSurface_Array_Test
        PATH MeshLib
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
        EXECUTABLE GocadTSurfaceReader
        EXECUTABLE_ARGS -i Top-Lower-Shaly.ts -o ${Data_BINARY_DIR}/MeshLib -b
        DEPENDS GocadTSurfaceReader-GocadTSurface_Mesh_Test-vtkdiff-mesh
        TESTER vtkdiff
        DIFF_DATA
        Top-Lower-Shaly.vtu Top-Lower-Shaly.vtu Reshape_Thickness Reshape_Thickness 1e-16 0
        Top-Lower-Shaly.vtu Top-Lower-Shaly.vtu Measured_Depth Measured_Depth 1e-16 0
    )
endif()

AddTest(
    NAME createIntermediateRasters_test
    PATH MeshGeoToolsLib/Hamburg
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Hamburg
    EXECUTABLE createIntermediateRasters
    EXECUTABLE_ARGS --file1 layer04.asc --file2 layer17.asc -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Hamburg/output.asc
    TESTER diff
    DIFF_DATA output0.asc
)

AddTest(
    NAME Vtu2Grid_Test
    PATH FileIO/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileIO
    EXECUTABLE Vtu2Grid
    EXECUTABLE_ARGS -i AmmerSubsurfaceCoarse.vtu -o ${Data_BINARY_DIR}/FileIO/AmmerGridOutput.vtu -x 200 -y 200 -z 20
    TESTER vtkdiff
    DIFF_DATA
    AmmerSubsurfaceGrid.vtu AmmerGridOutput.vtu MaterialIDs MaterialIDs 0 0
)

if(SNAKEMAKE AND TEE_TOOL_PATH)
    add_test(NAME snakemake_ExtractBoundary
        COMMAND ${SNAKEMAKE} -j 1
            --configfile ${PROJECT_BINARY_DIR}/buildinfo.yaml
            -s ${CMAKE_CURRENT_SOURCE_DIR}/ExtractBoundary.smk
    )

    add_test(NAME snakemake_VoxelGridFromLayers
        COMMAND ${SNAKEMAKE} -j 1
            --configfile ${PROJECT_BINARY_DIR}/buildinfo.yaml
            -s ${CMAKE_CURRENT_SOURCE_DIR}/VoxelGridFromLayers.smk
    )
    set_tests_properties(snakemake_ExtractBoundary snakemake_VoxelGridFromLayers
        PROPERTIES LABELS "default"
    )
    add_dependencies(ctest ExtractBoundary Layers2Grid AddFaultToVoxelGrid generateStructuredMesh)
endif()

AddTest(
    NAME partmesh_with_field_data
    PATH NodePartitionedMesh/partmesh
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/partmesh
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -n 2 -i cube_1x1x1_hex_8.vtu -x cube_1x1x1_hex_8 -o ${Data_BINARY_DIR}/NodePartitionedMesh/partmesh
    TESTER diff
    DIFF_DATA
        cube_1x1x1_hex_8_partitioned_cell_properties_cfg2.bin
        cube_1x1x1_hex_8_partitioned_cell_properties_val2.bin
        cube_1x1x1_hex_8_partitioned_msh_cfg2.bin
        cube_1x1x1_hex_8_partitioned_msh_ele2.bin
        cube_1x1x1_hex_8_partitioned_msh_ele_g2.bin
        cube_1x1x1_hex_8_partitioned_msh_nod2.bin
        cube_1x1x1_hex_8_partitioned_node_properties_cfg2.bin
        cube_1x1x1_hex_8_partitioned_node_properties_val2.bin
)

if(OGS_USE_NETCDF)
    AddTest(
        NAME NetCDF_2D_Test
        PATH FileConverter/
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileConverter
        EXECUTABLE NetCdfConverter
        EXECUTABLE_ARGS -i sresa1b_ncar_ccsm3-example.nc -o ${Data_BINARY_DIR}/FileConverter/sresa1b_ncar_ccsm3-example.vtu -v pr -t 0 --dim1 2 --dim2 1 --timestep-first 0 --timestep-last 0 -e tri
        TESTER vtkdiff
        DIFF_DATA
        sresa1b_ncar_ccsm3-example.vtu sresa1b_ncar_ccsm3-example.vtu pr pr 1e-16 0
    )

    AddTest(
        NAME NetCDF_3D_Test
        PATH FileConverter/
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileConverter
        EXECUTABLE NetCdfConverter
        EXECUTABLE_ARGS -i slim_100897_198.nc -o ${Data_BINARY_DIR}/FileConverter/slim_100897_198.vtu -v NO -t 0 --dim1 3 --dim2 2 --dim3 1 --timestep-first 0 --timestep-last 0 -e hex
        TESTER vtkdiff
        DIFF_DATA
        slim_100897_198.vtu slim_100897_198.vtu NO NO 1e-16 0
    )

    AddTest(
        NAME NetCDF_Image_Test
        PATH FileConverter
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileConverter
        EXECUTABLE NetCdfConverter
        EXECUTABLE_ARGS -i sresa1b_ncar_ccsm3-example.nc -o ${Data_BINARY_DIR}/FileConverter/sresa1b_ncar_ccsm3-example.asc -v pr -t 0 --dim1 2 --dim2 1 --timestep-first 0 --timestep-last 0 --images
        TESTER diff
        DIFF_DATA sresa1b_ncar_ccsm3-example0.asc
    )
endif()

AddTest(
    NAME RemoveGhostData_Test
    PATH MeshLib
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE RemoveGhostData
    EXECUTABLE_ARGS -i Mesh3D.pvtu -o ${Data_BINARY_DIR}/MeshLib/RemoveGhostDataOutput.vtu
    TESTER vtkdiff
    DIFF_DATA
    RemoveGhostDataOutput.vtu RemoveGhostDataOutput.vtu slice slice 0 0
)

AddTest(
    NAME RemoveGhostData_EllipticSquareTest
    PATH EllipticPETSc
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/EllipticPETSc
    EXECUTABLE RemoveGhostData
    EXECUTABLE_ARGS -i square_1e1_neumann_ts_1_t_1_000000.pvtu -o ${Data_BINARY_DIR}/EllipticPETSc/square_1e1_neumann_ts_1_t_1_000000.vtu
    TESTER diff
    DIFF_DATA
    square_1e1_neumann_ts_1_t_1_000000.vtu
)

AddTest(
    NAME Raster2Mesh_Elevation_Test
    PATH FileConverter
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileConverter
    EXECUTABLE Raster2Mesh
    EXECUTABLE_ARGS -i RainEvent30.asc -o ${Data_BINARY_DIR}/FileConverter/RainEvent30-elevation.vtu -e tri -p elevation
    TESTER diff
    DIFF_DATA RainEvent30-elevation.vtu
)

AddTest(
    NAME Raster2Mesh_Materials_Test
    PATH FileConverter
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileConverter
    EXECUTABLE Raster2Mesh
    EXECUTABLE_ARGS -i RainEvent30.asc -o ${Data_BINARY_DIR}/FileConverter/RainEvent30-materials.vtu -e quad -p materials
    TESTER vtkdiff
    DIFF_DATA
    RainEvent30-materials.vtu RainEvent30-materials.vtu MaterialIDs MaterialIDs 0 0
)

AddTest(
    NAME Raster2Mesh_Scalars_Test
    PATH FileConverter
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileConverter
    EXECUTABLE Raster2Mesh
    EXECUTABLE_ARGS -i RainEvent30.asc -o ${Data_BINARY_DIR}/FileConverter/RainEvent30-scalars.vtu -e tri -p scalar -n ScalarValues
    TESTER vtkdiff
    DIFF_DATA
    RainEvent30-scalars.vtu RainEvent30-scalars.vtu ScalarValues ScalarValues 0 0
)

AddTest(
    NAME AssignRasterDataToMesh2D_Test
    PATH MeshGeoToolsLib/Ammer
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Ammer
    EXECUTABLE AssignRasterDataToMesh
    EXECUTABLE_ARGS -i Ammer-Homogen100m-Final-TopSurface.vtu -r AmmerGWN.asc -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Ammer/AmmerGWN.vtu -s GWN -c -n
    TESTER vtkdiff
    DIFF_DATA
    AmmerGWN.vtu AmmerGWN.vtu GWN GWN 0 0
    AmmerGWN.vtu AmmerGWN.vtu GWN-2 GWN-2 0 0
)

AddTest(
    NAME AssignRasterDataToMesh1D_Test
    PATH MeshGeoToolsLib/Ammer
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Ammer
    EXECUTABLE AssignRasterDataToMesh
    EXECUTABLE_ARGS -i AmmerRivers.vtu -r AmmerGWN.asc -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Ammer/AmmerRiversGWN.vtu -s GWN -c -n
    TESTER vtkdiff
    DIFF_DATA
    AmmerRiversGWN.vtu AmmerRiversGWN.vtu GWN GWN 0 0
    AmmerRiversGWN.vtu AmmerRiversGWN.vtu GWN-2 GWN-2 0 0
)

AddTest(
    NAME ExtractMaterials_Test
    PATH MeshGeoToolsLib/Naegelstedt
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Naegelstedt
    EXECUTABLE ExtractMaterials
    EXECUTABLE_ARGS -i SmallTest.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Naegelstedt/SmallTest.vtu
    TESTER diff
    DIFF_DATA SmallTest_Layer1.vtu
              SmallTest_Layer2.vtu
              SmallTest_Layer3.vtu
)

# Tests requires gmsh
if(TARGET VerticalSliceFromLayers AND GMSH)
    AddTest(
        NAME VerticalSliceFromLayers_Test
        PATH MeshGeoToolsLib/Ammer
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Ammer
        EXECUTABLE VerticalSliceFromLayers
        EXECUTABLE_ARGS -i AmmerLayers.txt -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Ammer/AmmerSlice --start-x 3491000 --start-y 5380000 --end-x 3495000 --end-y 5385000 -r 200
        DIFF_DATA AmmerSlice.vtu AmmerSlice.vtu 1e-16
    )

endif()

if(TARGET GMSH2OGS AND SNAKEMAKE AND TEE_TOOL_PATH)
    add_test(NAME snakemake_GMSH2OGS_ExtractBoundary
        COMMAND ${SNAKEMAKE} --cores all
        --configfile ${PROJECT_BINARY_DIR}/buildinfo.yaml
        -s ${CMAKE_CURRENT_SOURCE_DIR}/GMSH2OGS_ExtractBoundary.smk
    )
    set_tests_properties(snakemake_GMSH2OGS_ExtractBoundary PROPERTIES LABELS "default")
    add_dependencies(ctest GMSH2OGS)
endif()

foreach(criterion ElementSize EdgeRatio EquiAngleSkew RadiusEdgeRatio SizeDifference)
    AddTest(
        NAME TrianglesGoodElementQuality_${criterion}_Test
        PATH MeshGeoToolsLib/Ammer
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Ammer
        EXECUTABLE AddElementQuality
        EXECUTABLE_ARGS -i AmmerGWN.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Ammer/AmmerGWNWithElementQuality_${criterion}.vtu -c ${criterion}
        TESTER vtkdiff
        DIFF_DATA
        AmmerGWNWithElementQuality.vtu AmmerGWNWithElementQuality_${criterion}.vtu ${criterion} ${criterion} 2e-7 2e-08
    )
endforeach()

foreach(criterion ElementSize EdgeRatio EquiAngleSkew RadiusEdgeRatio SizeDifference)
    AddTest(
        NAME TrianglesPoorElementQuality_${criterion}_Test
        PATH MeshGeoToolsLib/Hamburg
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Hamburg
        EXECUTABLE AddElementQuality
        EXECUTABLE_ARGS -i 00-surface.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Hamburg/00-surface-WithElementQuality_${criterion}.vtu -c ${criterion}
        TESTER vtkdiff
        DIFF_DATA
        00-surface-WithElementQuality.vtu 00-surface-WithElementQuality_${criterion}.vtu ${criterion} ${criterion} 2e-7 2e-8
    )
endforeach()

foreach(criterion ElementSize EdgeRatio EquiAngleSkew RadiusEdgeRatio SizeDifference)
    AddTest(
        NAME Mixed3DElementQuality_${criterion}_Test
        PATH FileIO
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileIO
        EXECUTABLE AddElementQuality
        EXECUTABLE_ARGS -i AmmerSubsurfaceCoarse.vtu -o ${Data_BINARY_DIR}/FileIO/AmmerSubsurfaceCoarse-WithElementQuality_${criterion}.vtu -c ${criterion}
        TESTER vtkdiff
        DIFF_DATA
        AmmerSubsurfaceCoarse-WithElementQuality.vtu AmmerSubsurfaceCoarse-WithElementQuality_${criterion}.vtu ${criterion} ${criterion} 1e-8 2e-11
    )
endforeach()

AddTest(
    NAME AddElementQuality_Behaelter_BE_ElementSize_Test
    PATH MeshGeoToolsLib/Behaelter_BE
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Behaelter_BE
    EXECUTABLE AddElementQuality
    EXECUTABLE_ARGS -i Behaelter_BE.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Behaelter_BE/Behaelter_BE_quality.vtu -c ElementSize
    TESTER vtkdiff
    DIFF_DATA
    Behaelter_BE_quality.vtu Behaelter_BE_quality.vtu ElementSize ElementSize 1e-12 1e-14
)

AddTest(
    NAME IntegrateBoreholesIntoMesh_MatOnly_Test
    PATH MeshGeoToolsLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib
    EXECUTABLE IntegrateBoreholesIntoMesh
    EXECUTABLE_ARGS -i PrismCube10x10x10.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/PrismBHE_mat.vtu -g testpoints.gml --min-id 4 --max-id 8
    TESTER diff
    DIFF_DATA
    PrismBHE_mat.vtu
)

AddTest(
    NAME IntegrateBoreholesIntoMesh_ElevationAndMat_Test
    PATH MeshGeoToolsLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib
    EXECUTABLE IntegrateBoreholesIntoMesh
    EXECUTABLE_ARGS -i PrismCube10x10x10.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/PrismBHE_elev.vtu -g testpoints.gml --min-id 4 --max-id 8 --min-elevation 4.5 --max-elevation 10
    TESTER diff
    DIFF_DATA
    PrismBHE_elev.vtu
)

AddTest(
    NAME ReviseMesh_Test
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE reviseMesh
    EXECUTABLE_ARGS -i basin_mesh.vtu -o ${Data_BINARY_DIR}/MeshLib/basin_mesh_fixed.vtu
    TESTER vtkdiff-mesh
    DIFF_DATA basin_mesh_fixed.vtu basin_mesh_fixed.vtu 1e-16
)

AddTest(
    NAME ReviseMesh_Test_Arrays
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE reviseMesh
    EXECUTABLE_ARGS -i basin_mesh.vtu -o ${Data_BINARY_DIR}/MeshLib/basin_mesh_fixed.vtu
    TESTER vtkdiff
    DIFF_DATA
    basin_mesh_fixed.vtu basin_mesh_fixed.vtu head head 0 0
    basin_mesh_fixed.vtu basin_mesh_fixed.vtu MaterialIDs MaterialIDs 0 0
)

if(TEST reviseMesh-ReviseMesh_Test-vtkdiff-mesh AND TEST reviseMesh-ReviseMesh_Test_Arrays)
    # Execute tests in order to prevent race condition
    set_tests_properties(reviseMesh-ReviseMesh_Test_Arrays PROPERTIES DEPENDS reviseMesh-ReviseMesh_Test-vtkdiff-mesh)
endif()

AddTest(
    NAME BinaryToPVTU
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/EllipticPETSc
    EXECUTABLE binaryToPVTU
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3 -o ${Data_BINARY_DIR}/EllipticPETSc/cube_1x1x1_hex_1e3
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex_1e3_0 cube_1x1x1_hex_1e3_0.vtu 1e-16
    cube_1x1x1_hex_1e3_1.vtu cube_1x1x1_hex_1e3_1.vtu 1e-16
    cube_1x1x1_hex_1e3_2.vtu cube_1x1x1_hex_1e3_2.vtu 1e-16
)

AddTest(
    NAME BinaryToPVTU_m1_load
    PATH NodePartitionedMesh/WithIntegrationPointStress
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress
    EXECUTABLE binaryToPVTU
    EXECUTABLE_ARGS -i m1_3Dload_ts_1_t_1_000000 -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/mi_load_binary_to_pvtu
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    mi_load_binary_to_pvtu_0.vtu mi_load_binary_to_pvtu_0.vtu 1e-16
    mi_load_binary_to_pvtu_1.vtu mi_load_binary_to_pvtu_1.vtu 1e-16
    mi_load_binary_to_pvtu_2.vtu mi_load_binary_to_pvtu_2.vtu 1e-16
)

AddTest(
    NAME BinaryToPVTU_pointheatsource
    PATH NodePartitionedMesh/WithIntegrationPointStress
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress
    EXECUTABLE binaryToPVTU
    EXECUTABLE_ARGS -i expected_pointheatsource_quadratic-mesh_ts_10_t_50000_000000 -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/pointheatsource_to_pvtu
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    pointheatsource_to_pvtu_0.vtu pointheatsource_to_pvtu_0.vtu 1e-16
    pointheatsource_to_pvtu_1.vtu pointheatsource_to_pvtu_1.vtu 1e-16
    pointheatsource_to_pvtu_2.vtu pointheatsource_to_pvtu_2.vtu 1e-16
    pointheatsource_to_pvtu_3.vtu pointheatsource_to_pvtu_3.vtu 1e-16
)

AddTest(
    NAME BinaryToPVTU_mesh_with_sigma_ip_and_mixed_element_types
    PATH NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    EXECUTABLE binaryToPVTU
    EXECUTABLE_ARGS -i mesh_with_3D_different_elements_sigma_ip -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/mesh_with_3D_different_elements_sigma_ip_pvtu
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    mesh_with_3D_different_elements_sigma_ip_pvtu_0.vtu mesh_with_3D_different_elements_sigma_ip_pvtu_0.vtu 1e-16
    mesh_with_3D_different_elements_sigma_ip_pvtu_1.vtu mesh_with_3D_different_elements_sigma_ip_pvtu_1.vtu 1e-16
)

AddTest(
    NAME BinaryToPVTU_mesh_with_sigma_ip_quad_tri
    PATH NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    EXECUTABLE binaryToPVTU
    EXECUTABLE_ARGS -i quad_tri_THM_t_864000_000000 -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad/quad_tri_THM_t_864000_000000_pvtu
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    quad_tri_THM_t_864000_000000_pvtu_0.vtu quad_tri_THM_t_864000_000000_pvtu_0.vtu 1e-16
    quad_tri_THM_t_864000_000000_pvtu_1.vtu quad_tri_THM_t_864000_000000_pvtu_1.vtu 1e-16
    quad_tri_THM_t_864000_000000_pvtu_2.vtu quad_tri_THM_t_864000_000000_pvtu_2.vtu 1e-16
    quad_tri_THM_t_864000_000000_pvtu_3.vtu quad_tri_THM_t_864000_000000_pvtu_3.vtu 1e-16
)

AddTest(
    NAME BinaryToPVTU_mesh_with_field_data_without_ip_data
    PATH NodePartitionedMesh/FieldDataWithoutIPData
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/FieldDataWithoutIPData
    EXECUTABLE binaryToPVTU
    EXECUTABLE_ARGS -i A2_tunnel_surface -o ${Data_BINARY_DIR}/NodePartitionedMesh/FieldDataWithoutIPData/A2_tunnel_surface_partitioned
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    A2_tunnel_surface_partitioned_0.vtu A2_tunnel_surface_partitioned_0.vtu 1e-16
    A2_tunnel_surface_partitioned_1.vtu A2_tunnel_surface_partitioned_1.vtu 1e-16
)

AddTest(
    NAME geometryToGmshAdaptiveGeo
    PATH MeshGeoToolsLib/geometryToGmshGeo/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/geometryToGmshGeo
    EXECUTABLE geometryToGmshGeo
    EXECUTABLE_ARGS -i square_1x1.gml -o ${Data_BINARY_DIR}/MeshGeoToolsLib/geometryToGmshGeo/square_1x1_adaptive.geo
    TESTER diff
    TESTER_ARGS --ignore-matching-lines=OpenGeoSys
    DIFF_DATA
    square_1x1_adaptive.geo
)

AddTest(
    NAME geometryToGmshHomogeneousGeo
    PATH MeshGeoToolsLib/geometryToGmshGeo/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/geometryToGmshGeo
    EXECUTABLE geometryToGmshGeo
    EXECUTABLE_ARGS -i square_1x1.gml -o ${Data_BINARY_DIR}/MeshGeoToolsLib/geometryToGmshGeo/square_1x1_homogeneous.geo --homogeneous -a 0.01
    TESTER diff
    TESTER_ARGS --ignore-matching-lines=OpenGeoSys
    DIFF_DATA
    square_1x1_homogeneous.geo
)

AddTest(
    NAME LineIntersectingDomainBoundary
    PATH MeshGeoToolsLib/geometryToGmshGeo/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/geometryToGmshGeo
    EXECUTABLE geometryToGmshGeo
    EXECUTABLE_ARGS -i square_1x1.gml -i line_intersecting_square.gml -o ${Data_BINARY_DIR}/MeshGeoToolsLib/geometryToGmshGeo/square_1x1_with_intersecting_line.geo
    PROPERTIES
        PASS_REGULAR_EXPRESSION
        "point with id 5 and coordinates \\(1.001000000001, 0.6, 0\\) is outside of the polygon"
)
AddTest(
    NAME ResetPropertiesInPolygonalRegion_AllElementNodesInPolygon
    PATH MeshGeoToolsLib/ResetPropertiesInPolygonalRegion/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/ResetPropertiesInPolygonalRegion
    EXECUTABLE ResetPropertiesInPolygonalRegion
    EXECUTABLE_ARGS -m Cube.vtu -n ValidCells -i 1 -g Polylines.gml -p Back -o ${Data_BINARY_DIR}/MeshGeoToolsLib/ResetPropertiesInPolygonalRegion/Cube-BackPolylinePropertyChange_all_element_nodes_inside.vtu
    TESTER vtkdiff
    DIFF_DATA
    Cube-BackPolylinePropertyChange_all_element_nodes_inside.vtu Cube-BackPolylinePropertyChange_all_element_nodes_inside.vtu ValidCells ValidCells 0 0
)

AddTest(
    NAME ResetPropertiesInPolygonalRegion_AtLeastOneElementNodeInPolygon
    PATH MeshGeoToolsLib/ResetPropertiesInPolygonalRegion/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/ResetPropertiesInPolygonalRegion
    EXECUTABLE ResetPropertiesInPolygonalRegion
    EXECUTABLE_ARGS -m Cube.vtu -n ValidCells -i 1 -g Polylines.gml -p Back --any_of -o ${Data_BINARY_DIR}/MeshGeoToolsLib/ResetPropertiesInPolygonalRegion/Cube-BackPolylinePropertyChange_at_least_one_element_node_inside.vtu
    TESTER vtkdiff
    DIFF_DATA
    Cube-BackPolylinePropertyChange_at_least_one_element_node_inside.vtu Cube-BackPolylinePropertyChange_at_least_one_element_node_inside.vtu ValidCells ValidCells 0 0
)

AddTest(
    NAME createRaster_10x10
    PATH GeoTools/createRaster/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/GeoTools/createRaster
    EXECUTABLE createRaster
    EXECUTABLE_ARGS -s 10 -r 10 -c 10 --ll_x 1000 --ll_y 100 -o ${Data_BINARY_DIR}/GeoTools/createRaster/10x10.asc
    TESTER diff
    DIFF_DATA
    10x10.asc 10x10.asc
)

AddTest(
    NAME addDataToRaster_10x10
    PATH GeoTools/addDataToRaster/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/GeoTools/addDataToRaster
    EXECUTABLE addDataToRaster
    EXECUTABLE_ARGS --ll_x 1000 --ll_y 100 --ur_x 1100 --ur_y 200 --function sinxsiny --scaling_value 1 --offset_value 0 -i ${Data_SOURCE_DIR}/GeoTools/createRaster/10x10.asc -o ${Data_BINARY_DIR}/GeoTools/addDataToRaster/10x10_sinxsiny.asc
    TESTER diff
    DIFF_DATA
    10x10_sinxsiny.asc 10x10_sinxsiny.asc
)

AddTest(
    NAME addDataToRaster_10x10_step_function
    PATH GeoTools/addDataToRaster/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/GeoTools/addDataToRaster
    EXECUTABLE addDataToRaster
    EXECUTABLE_ARGS --ll_x 1000 --ll_y 100 --ur_x 1100 --ur_y 200 --function step --scaling_value 1 --offset_value 0 -i ${Data_SOURCE_DIR}/GeoTools/createRaster/10x10.asc -o ${Data_BINARY_DIR}/GeoTools/addDataToRaster/10x10_step.asc
    TESTER diff
    DIFF_DATA
    10x10_step.asc 10x10_step.asc
)

AddTest(
    NAME GMSH2OGS_linearElements
    PATH Utils/GMSH2OGS
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Utils/GMSH2OGS
    EXECUTABLE GMSH2OGS
    EXECUTABLE_ARGS -i linear_mesh.msh
                    -o ${Data_BINARY_DIR}/Utils/GMSH2OGS/linear_mesh.vtu
    DIFF_DATA linear_mesh.vtu linear_mesh.vtu 1.e-16
)

AddTest(
    NAME GMSH2OGS_quadratic_quadrilateral
    PATH Utils/GMSH2OGS
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Utils/GMSH2OGS
    EXECUTABLE GMSH2OGS
    EXECUTABLE_ARGS -i quadratic_quadrilateral.msh
                    -o ${Data_BINARY_DIR}/Utils/GMSH2OGS/quadratic_quadrilateral.vtu
    TESTER vtkdiff-mesh
    DIFF_DATA quadratic_quadrilateral.vtu quadratic_quadrilateral.vtu 1.e-16
)

AddTest(
    NAME GMSH2OGS_quadratic_elements
    PATH Utils/GMSH2OGS
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Utils/GMSH2OGS
    EXECUTABLE GMSH2OGS
    EXECUTABLE_ARGS -i quadratic_mesh.msh
                    -o ${Data_BINARY_DIR}/Utils/GMSH2OGS/quadratic_mesh.vtu
    REQUIREMENTS NOT (OGS_USE_MPI)
    TESTER vtkdiff-mesh
    DIFF_DATA quadratic_mesh.vtu quadratic_mesh.vtu 1.e-16
)

AddTest(
    NAME generateGeometry_10x10_quad
    PATH GeoTools/generateGeometry/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/GeoTools/generateGeometry
    EXECUTABLE generateGeometry
    EXECUTABLE_ARGS --x0 0 --y0 0 --z0 0 --x1 10 --y1 10 --z1 0 --nx 1 --ny 1 --geometry_name TestGeometry --polyline_name TestQuad -o ${Data_BINARY_DIR}/GeoTools/generateGeometry/TestGeometry_10x10_quad.gml
    TESTER diff
    DIFF_DATA
    TestGeometry_10x10_quad.gml TestGeometry_10x10_quad.gml
)

AddTest(
    NAME generateGeometry_10_line_nx9
    PATH GeoTools/generateGeometry/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/GeoTools/generateGeometry
    EXECUTABLE generateGeometry
    EXECUTABLE_ARGS --x0 0 --y0 0 --z0 0 --x1 10 --y1 0 --z1 0 --nx 9 --geometry_name TestGeometry --polyline_name TestLine -o ${Data_BINARY_DIR}/GeoTools/generateGeometry/TestGeometry_10_line_nx9.gml
    TESTER diff
    DIFF_DATA
    TestGeometry_10_line_nx9.gml TestGeometry_10_line_nx9.gml
)

AddTest(
    NAME generateGeometry_point
    PATH GeoTools/generateGeometry/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/GeoTools/generateGeometry
    EXECUTABLE generateGeometry
    EXECUTABLE_ARGS --x0 1 --y0 1 --z0 0 --x1 1 --y1 1 --z1 0 --geometry_name TestGeometry --polyline_name TestPoint  -o ${Data_BINARY_DIR}/GeoTools/generateGeometry/TestGeometry_point.gml
    TESTER diff
    DIFF_DATA
    TestGeometry_point.gml TestGeometry_point.gml
)

AddTest(
    NAME NodeOrdering_M0
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE NodeReordering
    EXECUTABLE_ARGS -i ReorderTestMesh.vtu -o ${Data_BINARY_DIR}/MeshLib/ReorderTestMeshM0.vtu -m 0
    REQUIREMENTS NOT (OGS_USE_MPI)
    TESTER vtkdiff-mesh
    DIFF_DATA ReorderTestMeshM0.vtu ReorderTestMeshM0.vtu 1.e-16
)

AddTest(
    NAME NodeOrdering_M1
    PATH MeshLib/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE NodeReordering
    EXECUTABLE_ARGS -i ReorderTestMesh.vtu -o ${Data_BINARY_DIR}/MeshLib/ReorderTestMeshM1.vtu -m 1
    REQUIREMENTS NOT (OGS_USE_MPI)
    TESTER vtkdiff-mesh
    DIFF_DATA ReorderTestMeshM1.vtu ReorderTestMeshM1.vtu 1.e-16
)

AddTest(
    NAME 1D_HeatConduction_dirichlet-line_60_heat
    PATH Parabolic/T/1D_dirichlet
    EXECUTABLE PVD2XDMF
    EXECUTABLE_ARGS ${Data_SOURCE_DIR}/Parabolic/T/1D_dirichlet/line_60_heat.pvd
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER xdmfdiff
    DIFF_DATA
    line_60_heat_line_60_heat_ts_0_t_0.000000.xdmf line_60_heat_line_60_heat_ts_0_t_0.000000.xdmf MaterialIDs MaterialIDs 0 0
    line_60_heat_line_60_heat_ts_0_t_0.000000.xdmf line_60_heat_line_60_heat_ts_0_t_0.000000.xdmf HeatFlowRate HeatFlowRate 0 0
    line_60_heat_line_60_heat_ts_0_t_0.000000.xdmf line_60_heat_line_60_heat_ts_0_t_0.000000.xdmf temperature temperature 1e-13 0
    line_60_heat_line_60_heat_ts_0_t_0.000000.xdmf line_60_heat_line_60_heat_ts_0_t_0.000000.xdmf heat_flux heat_flux 1e-13 0
)

AddTest(
    NAME MapMeshToMesh_Test
    PATH MeshGeoToolsLib/Hamburg
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Hamburg
    EXECUTABLE MeshMapping
    EXECUTABLE_ARGS -i plain.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Hamburg/meshmapping.vtu -m 00-surface.vtu -d 150
    TESTER vtkdiff-mesh
    DIFF_DATA meshmapping.vtu meshmapping.vtu 8e-3
)

AddTest(
    NAME MapMeshToRasterASC_Test
    PATH MeshGeoToolsLib/Hamburg
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshGeoToolsLib/Hamburg
    EXECUTABLE MeshMapping
    EXECUTABLE_ARGS -i plain.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Hamburg/rastermapping.vtu -r 00-raster.asc -s 100
    TESTER vtkdiff-mesh
    DIFF_DATA rastermapping.vtu rastermapping.vtu 1.5e-14
)

AddTest(
    NAME MapMeshToRasterXYZ_Test
    PATH FileIO
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileIO
    EXECUTABLE MeshMapping
    EXECUTABLE_ARGS -i XyzTest.vtu -o ${Data_BINARY_DIR}/FileIO/XyzTest-Mapped.vtu -r XyzTest.xyz
    TESTER vtkdiff-mesh
    DIFF_DATA XyzTest-Mapped.vtu XyzTest-Mapped.vtu 1e-12
)

AddTest(
    NAME Utils_pvtu2vtu
    PATH Utils/PVTU2VTU
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Utils/PVTU2VTU
    EXECUTABLE pvtu2vtu
    EXECUTABLE_ARGS -i quad_tri_THM_MPI_t_86400_000000.pvtu -o ${Data_BINARY_DIR}/Utils/PVTU2VTU/merged_quad_tri_THM_MPI_t_86400_000000.vtu
    TESTER vtkdiff
    DIFF_DATA
    merged_quad_tri_THM_MPI_t_86400_000000.vtu merged_quad_tri_THM_MPI_t_86400_000000.vtu displacement displacement 1e-15 1e-15
    merged_quad_tri_THM_MPI_t_86400_000000.vtu merged_quad_tri_THM_MPI_t_86400_000000.vtu pressure pressure 1e-15 1e-15
    merged_quad_tri_THM_MPI_t_86400_000000.vtu merged_quad_tri_THM_MPI_t_86400_000000.vtu temperature temperature 1e-15 1e-15
    merged_quad_tri_THM_MPI_t_86400_000000.vtu merged_quad_tri_THM_MPI_t_86400_000000.vtu epsilon epsilon 1e-15 1e-15
    merged_quad_tri_THM_MPI_t_86400_000000.vtu merged_quad_tri_THM_MPI_t_86400_000000.vtu sigma sigma 1e-15 1e-15
)

AddTest(
    NAME mixedElements
    PATH NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements
    EXECUTABLE ipDataToPointCloud
    EXECUTABLE_ARGS -i mesh_with_3D_different_elements_sigma_ip.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/mesh_with_3D_different_elements_sigma_ip_point_cloud.vtu
    TESTER vtkdiff
    DIFF_DATA
    mesh_with_3D_different_elements_sigma_ip_point_cloud.vtu mesh_with_3D_different_elements_sigma_ip_point_cloud.vtu sigma_ip sigma_ip 1e-15 0
)

AddTest(
    NAME triAndQuadMesh
    PATH NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad
    EXECUTABLE ipDataToPointCloud
    EXECUTABLE_ARGS -i quad_tri_THM_t_864000_000000.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/TriQuad/quad_tri_THM_t_864000_000000_point_cloud.vtu
    TESTER vtkdiff
    DIFF_DATA
    quad_tri_THM_t_864000_000000_point_cloud.vtu quad_tri_THM_t_864000_000000_point_cloud.vtu sigma_ip sigma_ip 1e-15 1e-15
    quad_tri_THM_t_864000_000000_point_cloud.vtu quad_tri_THM_t_864000_000000_point_cloud.vtu epsilon_ip epsilon_ip 1e-15 0
)

AddTest(
    NAME m1_3Dsquare
    PATH Mechanics/m1_3Dsquare
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Mechanics/m1_3Dsquare
    EXECUTABLE ipDataToPointCloud
    EXECUTABLE_ARGS -i m1_3Dsquare_ts_1_t_1.000000.vtu -o ${Data_BINARY_DIR}/Mechanics/m1_3Dsquare/m1_3Dsquare_ts_1_t_1.000000_point_cloud.vtu
    TESTER vtkdiff
    DIFF_DATA
    m1_3Dsquare_ts_1_t_1.000000_point_cloud.vtu m1_3Dsquare_ts_1_t_1.000000_point_cloud.vtu sigma_ip sigma_ip 1e-15 0
)

if(OGS_BUILD_SWMM)
    AddTest(
        NAME SWMM_INP_Geo_Test
        PATH FileConverter
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileConverter
        EXECUTABLE SWMMConverter
        EXECUTABLE_ARGS -i TestExample_SC2.inp -g ${Data_BINARY_DIR}/FileConverter/TestExample_SC2.gml
        TESTER diff
        DIFF_DATA
        TestExample_SC2.gml TestExample_SC2.gml
    )

    AddTest(
        NAME SWMM_INP_Mesh_Test
        PATH FileConverter
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileConverter
        EXECUTABLE SWMMConverter
        EXECUTABLE_ARGS -i TestExample_SC2.inp -m ${Data_BINARY_DIR}/FileConverter/TestExample_SC2.vtu
        TESTER vtkdiff-mesh
        DIFF_DATA TestExample_SC2.vtu TestExample_SC2.vtu 1.e-16
    )
endif()

AddTest(
    NAME xyz2asc_Test
    PATH FileIO/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileIO
    EXECUTABLE Raster2ASC
    EXECUTABLE_ARGS -i XyzTest.xyz -o ${Data_BINARY_DIR}/FileIO/XyzTest.asc
    TESTER diff
    DIFF_DATA XyzTest.asc
)

AddTest(
    NAME grd2asc_Test
    PATH FileIO/
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileIO
    EXECUTABLE Raster2ASC
    EXECUTABLE_ARGS -i GrdTest.grd -o ${Data_BINARY_DIR}/FileIO/GrdTest.asc
    TESTER diff
    DIFF_DATA GrdTest.asc
)

AddTest(
    NAME RemoveUnusedPoints_Cube
    PATH FileIO/RemoveUnusedPoints
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileIO/RemoveUnusedPoints
    EXECUTABLE RemoveUnusedPoints
    EXECUTABLE_ARGS -i cube_1x1x1_with_additional_points.gml -o ${Data_BINARY_DIR}/FileIO/RemoveUnusedPoints/cube_1x1x1_with_additional_points_cleaned.gml
    TESTER diff
    DIFF_DATA cube_1x1x1_with_additional_points_cleaned.gml cube_1x1x1_with_additional_points_cleaned.gml
)

AddTest(
    NAME RemoveUnusedPoints_WESSRivers
    PATH FileIO/RemoveUnusedPoints
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/FileIO/RemoveUnusedPoints
    EXECUTABLE RemoveUnusedPoints
    EXECUTABLE_ARGS -i WESSRivers.gml -o ${Data_BINARY_DIR}/FileIO/RemoveUnusedPoints/WESSRivers_cleaned.gml
    TESTER diff
    DIFF_DATA WESSRivers_cleaned.gml WESSRivers_cleaned.gml
)

AddTest(
    NAME ConstructMeshesFromGeometry_square_lines_and_points
    PATH MeshGeoToolsLib/ConstructMeshesFromGeometry
    WORKING_DIRECTORY ${Data_BINARY_DIR}/<PATH>
    EXECUTABLE constructMeshesFromGeometry
    EXECUTABLE_ARGS -m ${Data_SOURCE_DIR}/<PATH>/square_1x1_quad8_1e2.vtu -g ${Data_SOURCE_DIR}/<PATH>/square_1x1.gml
    TESTER vtkdiff-mesh
    DIFF_DATA
    square_1x1_geometry_origin.vtu square_1x1_geometry_origin.vtu 1e-16
    square_1x1_geometry_left.vtu   square_1x1_geometry_left.vtu 1e-16
    square_1x1_geometry_right.vtu  square_1x1_geometry_right.vtu 1e-16
    square_1x1_geometry_bottom.vtu square_1x1_geometry_bottom.vtu 1e-16
    square_1x1_geometry_top.vtu    square_1x1_geometry_top.vtu 1e-16
)

if(NOT OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE ../../web/content/docs/tutorials/bhe_meshing/notebook-bhe_meshing.md
                 PYTHON_PACKAGES openpyxl
                 RUNTIME 10)
    NotebookTest(NOTEBOOKFILE ../../web/content/docs/tutorials/mesh_tutorial/notebook-mesh_tutorial.md
                 RUNTIME 10)
    NotebookTest(NOTEBOOKFILE ../../web/content/docs/tutorials/Inclined_bhe_meshing/notebook-inclined_bhe_meshing.md
                 RUNTIME 10)
endif()

AddTest(
    NAME RemoveMeshElements_AABB_2D_regular
    PATH MeshLib
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE removeMeshElements
    EXECUTABLE_ARGS -i AREHS_Layer17.vtu
                    -o ${Data_BINARY_DIR}/MeshLib/AREHS_2D_AABB_regular.vtu
                    --x-min 12000 --x-max 15000 --y-min 12000
    REQUIREMENTS NOT (OGS_USE_MPI)
    TESTER vtkdiff-mesh
    DIFF_DATA AREHS_2D_AABB_regular.vtu AREHS_2D_AABB_regular.vtu 1.e-16
)

AddTest(
    NAME RemoveMeshElements_AABB_2D_inverted
    PATH MeshLib
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE removeMeshElements
    EXECUTABLE_ARGS -i AREHS_Layer17.vtu
                    -o ${Data_BINARY_DIR}/MeshLib/AREHS_2D_AABB_inverted.vtu
                    --x-min 12000 --x-max 15000 --y-min 12000 --invert
    REQUIREMENTS NOT (OGS_USE_MPI)
    TESTER vtkdiff-mesh
    DIFF_DATA AREHS_2D_AABB_inverted.vtu AREHS_2D_AABB_inverted.vtu 1.e-16
)

AddTest(
    NAME RemoveMeshElements_AABB_3D_regular
    PATH MeshLib
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE removeMeshElements
    EXECUTABLE_ARGS -i AREHS_test.vtu
                    -o ${Data_BINARY_DIR}/MeshLib/AREHS_3D_AABB_regular.vtu
                    --x-min 12000 --x-max 15000 --y-min 12000 --z-min -3000 --z-max -2000
    REQUIREMENTS NOT (OGS_USE_MPI)
    TESTER vtkdiff-mesh
    DIFF_DATA AREHS_3D_AABB_regular.vtu AREHS_3D_AABB_regular.vtu 1.e-16
)

AddTest(
    NAME RemoveMeshElements_AABB_3D_inverted
    PATH MeshLib
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/MeshLib
    EXECUTABLE removeMeshElements
    EXECUTABLE_ARGS -i AREHS_test.vtu
                    -o ${Data_BINARY_DIR}/MeshLib/AREHS_3D_AABB_inverted.vtu
                    --x-min 12000 --x-max 15000 --y-min 12000 --z-min -3000 --z-max -2000 --invert
    REQUIREMENTS NOT (OGS_USE_MPI)
    TESTER vtkdiff-mesh
    DIFF_DATA AREHS_3D_AABB_inverted.vtu AREHS_3D_AABB_inverted.vtu 1.e-16
)

if(OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE Utils/partmesh/partmesh_roundtrip.md RUNTIME 10 SKIP_WEB)
endif()
