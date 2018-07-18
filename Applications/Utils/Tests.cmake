
AddTest(
    NAME MapGeometryToMeshSurface_Ammer
    PATH MeshGeoToolsLib/Ammer
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m Ammer-Homogen100m-Final-TopSurface.vtu -i Ammer-Rivers.gml -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Ammer/Ammer-Rivers-Mapped.gml
    TESTER diff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA Ammer-Rivers-Mapped.gml
)

AddTest(
    NAME MapGeometryToMeshSurface_Bode
    PATH MeshGeoToolsLib/Bode
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m BodeComplex.msh -i BodeEZG_Fliessgewaesser.gml -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Bode/BodeEZG_Fliessgewaesser-Mapped.gml
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER diff
    DIFF_DATA BodeEZG_Fliessgewaesser-Mapped.gml
)

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
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    DIFF_DATA
    expected_post_single_joint_pcs_0_ts_1_t_1.000000.vtu post_single_joint_pcs_0_ts_1_t_1.000000.vtu u u 1e-14 1e-14
)

# Mac is producing slightly different partitioning, so the results are not
# comparable.
AddTest(
    NAME partmesh_2Dmesh_3partitions_ascii
    PATH NodePartitionedMesh/partmesh_2Dmesh_3partitions
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -a -m -n 3 -i 2Dmesh.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/partmesh_2Dmesh_3partitions
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
    PATH NodePartitionedMesh/partmesh_2Dmesh_3partitions
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -m -n 3 -i 2Dmesh.vtu
                    -o ${Data_BINARY_DIR}/NodePartitionedMesh/partmesh_2Dmesh_3partitions --
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
)
