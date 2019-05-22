
AddTest(
    NAME MapGeometryToMeshSurface_Ammer
    PATH MeshGeoToolsLib/Ammer
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m Ammer-Homogen100m-Final-TopSurface.vtu -i Ammer-Rivers.gml -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Ammer/Ammer-Rivers-Mapped.gml
    TESTER diff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA Ammer-Rivers-Mapped.gml
)

# Disable test on eve frontends
if(NOT "${HOSTNAME}" MATCHES "frontend.*")
    AddTest(
        NAME MapGeometryToMeshSurface_Bode
        PATH MeshGeoToolsLib/Bode
        EXECUTABLE MapGeometryToMeshSurface
        EXECUTABLE_ARGS -m BodeComplex.msh -i BodeEZG_Fliessgewaesser.gml -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Bode/BodeEZG_Fliessgewaesser-Mapped.gml
        REQUIREMENTS NOT OGS_USE_MPI
        TESTER diff
        DIFF_DATA BodeEZG_Fliessgewaesser-Mapped.gml
    )
endif()

AddTest(
    NAME MapGeometryToMeshSurface_Naegelstedt
    PATH MeshGeoToolsLib/Naegelstedt
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m SmallTest.vtu -i RiverNetwork.gml -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Naegelstedt/RiverNetwork-Mapped.gml
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER diff
    DIFF_DATA RiverNetwork-Mapped.gml
)

AddTest(
    NAME postLIE
    PATH LIE/PostProcessing
    EXECUTABLE postLIE
    EXECUTABLE_ARGS -i single_joint_pcs_0.pvd -o ${Data_BINARY_DIR}/LIE/PostProcessing/post_single_joint_pcs_0.pvd
    REQUIREMENTS NOT OGS_USE_MPI AND OGS_BUILD_PROCESS_LIE
    TESTER vtkdiff
    DIFF_DATA
    expected_post_single_joint_pcs_0_ts_1_t_1.000000.vtu post_single_joint_pcs_0_ts_1_t_1.000000.vtu u u 1e-14 1e-14
)

AddTest(
    NAME identifySubdomains_2D_Create
    PATH MeshGeoToolsLib/IdentifySubdomains
    EXECUTABLE identifySubdomains
    EXECUTABLE_ARGS -m 2D_mesh.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/IdentifySubdomains/new_ -- 2D_mesh_top_boundary.vtu 2D_mesh_bottom_boundary.vtu
    REQUIREMENTS NOT OGS_USE_MPI
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
    EXECUTABLE identifySubdomains
    EXECUTABLE_ARGS -m 2D_mesh.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/IdentifySubdomains/check_ -- 2D_mesh_top.vtu 2D_mesh_bottom.vtu
    REQUIREMENTS NOT OGS_USE_MPI
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
    EXECUTABLE identifySubdomains
    EXECUTABLE_ARGS -m river_domain_triangle.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/IdentifySubdomains/triangle_ -- river_bc.vtu
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    DIFF_DATA
    river_bc_triangle.vtu triangle_river_bc.vtu bulk_node_ids bulk_node_ids 0 0
    #river_bc_triangle.vtu triangle_river_bc.vtu bulk_element_ids bulk_element_ids 0 0   # TODO (naumov) Needs extension of vtkdiff to FieldData
    river_bc_triangle.vtu triangle_river_bc.vtu number_bulk_elements number_bulk_elements 0 0
)

AddTest(
    NAME identifySubdomains_riverPrismMesh
    PATH MeshGeoToolsLib/IdentifySubdomains
    EXECUTABLE identifySubdomains
    EXECUTABLE_ARGS -s 1e-3 -m river_domain_prism.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/IdentifySubdomains/prism_ -- river_bc.vtu
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    DIFF_DATA
    river_bc_prism.vtu prism_river_bc.vtu bulk_node_ids bulk_node_ids 0 0
    #river_bc_prism.vtu prism_river_bc.vtu bulk_element_ids bulk_element_ids 0 0 # TODO (naumov) Needs extension of vtkdiff to FieldData
    river_bc_prism.vtu prism_river_bc.vtu number_bulk_elements number_bulk_elements 0 0
)

# Mac is producing slightly different partitioning, so the results are not
# comparable.
AddTest(
    NAME partmesh_2Dmesh_3partitions_ascii
    PATH NodePartitionedMesh/partmesh_2Dmesh_3partitions/ASCII
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -a -m -n 3 -i 2Dmesh.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/partmesh_2Dmesh_3partitions/ASCII
    REQUIREMENTS NOT (OGS_USE_MPI OR APPLE)
    TESTER diff
    DIFF_DATA 2Dmesh_partitioned_elems_3.msh
              2Dmesh_partitioned_cfg3.msh
              2Dmesh_partitioned_nodes_3.msh
)

# Mac is producing slightly different partitioning, so the results are not
# comparable.
AddTest(
    NAME partmesh_2Dmesh_3partitions_binary
    PATH NodePartitionedMesh/partmesh_2Dmesh_3partitions/Binary
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 3 -i 2Dmesh.vtu
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/partmesh_2Dmesh_3partitions/Binary --
                    2Dmesh_PLY_EAST.vtu
                    2Dmesh_PLY_WEST.vtu
                    2Dmesh_PLY_NORTH.vtu
                    2Dmesh_PLY_SOUTH.vtu
                    2Dmesh_POINT4.vtu
                    2Dmesh_POINT5.vtu
    REQUIREMENTS NOT (OGS_USE_MPI OR APPLE)
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

# Regression test for https://github.com/ufz/ogs/issues/1845 fixed in
# https://github.com/ufz/ogs/pull/2237
# checkMesh crashed when encountered Line3 element.
AddTest(
    NAME checkMesh_LIE_HM_TaskB
    PATH LIE/HydroMechanics
    EXECUTABLE checkMesh
    EXECUTABLE_ARGS -p -v TaskB_mesh.vtu
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME mesh2raster_test
    PATH MeshGeoToolsLib/Hamburg
    EXECUTABLE Mesh2Raster
    EXECUTABLE_ARGS -i 00-surface.vtu -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Hamburg/00-raster.asc -c 25
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER diff
    DIFF_DATA 00-raster.asc
)

MeshTest(
    NAME ExtractSurfaceLeft
    PATH MeshLib/
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3_layers_10.vtu -o ${Data_BINARY_DIR}/MeshLib/Left.vtu -x 1 -y 0 -z 0 -a 25
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA Left.vtu Left.vtu 1e-16
)

MeshTest(
    NAME ExtractSurfaceRight
    PATH MeshLib/
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3_layers_10.vtu -o ${Data_BINARY_DIR}/MeshLib/Right.vtu -x -1 -y 0 -z 0 -a 25
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA Right.vtu Right.vtu 1e-16
)

MeshTest(
    NAME ExtractSurfaceFront
    PATH MeshLib/
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3_layers_10.vtu -o ${Data_BINARY_DIR}/MeshLib/Front.vtu -x 0 -y 1 -z 0 -a 25
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA Front.vtu Front.vtu 1e-16
)

MeshTest(
    NAME ExtractSurfaceBack
    PATH MeshLib/
    EXECUTABLE ExtractSurface
    EXECUTABLE_ARGS -i cube_1x1x1_hex_1e3_layers_10.vtu -o ${Data_BINARY_DIR}/MeshLib/Back.vtu -x 0 -y -1 -z 0 -a 25
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA Back.vtu Back.vtu 1e-16
)
