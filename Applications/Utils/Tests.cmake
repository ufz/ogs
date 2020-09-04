
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

MeshTest(
    NAME GocadTSurface_Mesh_Test
    PATH MeshLib/
    EXECUTABLE GocadTSurfaceReader
    EXECUTABLE_ARGS -i Top-Lower-Shaly.ts -o ${Data_BINARY_DIR}/MeshLib -b
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA Top-Lower-Shaly.vtu Top-Lower-Shaly.vtu 1e-16
)

AddTest(
    NAME GocadTSurface_Array_Test
    PATH MeshLib/
    EXECUTABLE GocadTSurfaceReader
    EXECUTABLE_ARGS -i Top-Lower-Shaly.ts -o ${Data_BINARY_DIR}/MeshLib -b
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    DIFF_DATA
    Top-Lower-Shaly.vtu Top-Lower-Shaly.vtu Reshape_Thickness Reshape_Thickness 1e-16 0
    Top-Lower-Shaly.vtu Top-Lower-Shaly.vtu Measured_Depth Measured_Depth 1e-16 0
)

AddTest(
    NAME createIntermediateRasters_test
    PATH MeshGeoToolsLib/Hamburg
    EXECUTABLE createIntermediateRasters
    EXECUTABLE_ARGS --file1 layer04.asc --file2 layer17.asc -o ${Data_BINARY_DIR}/MeshGeoToolsLib/Hamburg/output.asc
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER diff
    DIFF_DATA output0.asc
)

AddTest(
    NAME Vtu2Grid_Test
    PATH FileIO/
    EXECUTABLE Vtu2Grid
    EXECUTABLE_ARGS -i AmmerSubsurfaceCoarse.vtu -o ${Data_BINARY_DIR}/FileIO/AmmerGridOutput.vtu -x 200 -y 200 -z 20
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    DIFF_DATA
    AmmerSubsurfaceGrid.vtu AmmerGridOutput.vtu MaterialIDs MaterialIDs 0 0
)

foreach(element_type tri quad)
    if(OGS_USE_MPI)
        break()
    endif()
    add_test(NAME ExtractBoundary-generate_input_files_${element_type}
        COMMAND generateStructuredMesh -e ${element_type} --lx 1 --ly 1 --nx 10 --ny 10 -o ${Data_BINARY_DIR}/FileIO/square_1x1_${element_type}.vtu
        WORKING_DIRECTORY ${Data_BINARY_DIR}/FileIO/
    )
    set_tests_properties(ExtractBoundary-generate_input_files_${element_type}
        PROPERTIES COST 1) # RUNTIME 1

    AddTest(
        NAME ExtractBoundary_2D_${element_type}_Test
        PATH FileIO/
        EXECUTABLE ExtractBoundary
        EXECUTABLE_ARGS -i ${Data_BINARY_DIR}/FileIO/square_1x1_${element_type}.vtu -o ${Data_BINARY_DIR}/FileIO/square_1x1_${element_type}_boundary.vtu
        REQUIREMENTS NOT OGS_USE_MPI
        TESTER vtkdiff
        DIFF_DATA
        square_1x1_${element_type}_boundary.vtu square_1x1_${element_type}_boundary.vtu bulk_node_ids bulk_node_ids 0 0
        square_1x1_${element_type}_boundary.vtu square_1x1_${element_type}_boundary.vtu bulk_element_ids bulk_element_ids 0 0
        square_1x1_${element_type}_boundary.vtu square_1x1_${element_type}_boundary.vtu bulk_face_ids bulk_face_ids 0 0
        RUNTIME 1
        DEPENDS ExtractBoundary-generate_input_files_${element_type}
    )
endforeach()

AddTest(
    NAME partmesh_with_field_data
    PATH NodePartitionedMesh/partmesh
    EXECUTABLE partmesh
    EXECUTABLE_ARGS -n 2 -i cube_1x1x1_hex_8.vtu -o ${Data_BINARY_DIR}/NodePartitionedMesh/partmesh
    TESTER diff
    REQUIREMENTS NOT OGS_USE_MPI
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
        EXECUTABLE NetCdfConverter
        EXECUTABLE_ARGS -i sresa1b_ncar_ccsm3-example.nc -o ${Data_BINARY_DIR}/FileConverter/sresa1b_ncar_ccsm3-example.vtu -v pr -t 0 --dim1 2 --dim2 1 --timestep-first 0 --timestep-last 0 -e tri
        REQUIREMENTS NOT OGS_USE_MPI
        TESTER vtkdiff
        DIFF_DATA
        sresa1b_ncar_ccsm3-example.vtu sresa1b_ncar_ccsm3-example.vtu pr pr 1e-16 0
    )

    AddTest(
        NAME NetCDF_3D_Test
        PATH FileConverter/
        EXECUTABLE NetCdfConverter
        EXECUTABLE_ARGS -i slim_100897_198.nc -o ${Data_BINARY_DIR}/FileConverter/slim_100897_198.vtu -v NO -t 0 --dim1 3 --dim2 2 --dim3 1 --timestep-first 0 --timestep-last 0 -e hex
        REQUIREMENTS NOT OGS_USE_MPI
        TESTER vtkdiff
        DIFF_DATA
        slim_100897_198.vtu slim_100897_198.vtu NO NO 1e-16 0
    )

    AddTest(
        NAME NetCDF_Image_Test
        PATH FileConverter
        EXECUTABLE NetCdfConverter
        EXECUTABLE_ARGS -i sresa1b_ncar_ccsm3-example.nc -o ${Data_BINARY_DIR}/FileConverter/sresa1b_ncar_ccsm3-example.asc -v pr -t 0 --dim1 2 --dim2 1 --timestep-first 0 --timestep-last 0 --images
        REQUIREMENTS NOT OGS_USE_MPI
        TESTER diff
        DIFF_DATA sresa1b_ncar_ccsm3-example0.asc
    )
endif()

if(OGS_BUILD_GUI)
    AddTest(
        NAME RemoveGhostData_Test
        PATH MeshLib/
        EXECUTABLE RemoveGhostData
        EXECUTABLE_ARGS -i Mesh3D.pvtu -o ${Data_BINARY_DIR}/MeshLib/RemoveGhostDataOutput.vtu
        REQUIREMENTS NOT OGS_USE_MPI
        TESTER vtkdiff
        DIFF_DATA
        RemoveGhostDataOutput.vtu RemoveGhostDataOutput.vtu slice slice 0 0
    )
endif()

AddTest(
    NAME Raster2Mesh_Elevation_Test
    PATH FileConverter
    EXECUTABLE Raster2Mesh
    EXECUTABLE_ARGS -i RainEvent30.asc -o ${Data_BINARY_DIR}/FileConverter/RainEvent30-elevation.vtu -e tri -p elevation
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER diff
    DIFF_DATA RainEvent30-elevation.vtu
)

AddTest(
    NAME Raster2Mesh_Materials_Test
    PATH FileConverter
    EXECUTABLE Raster2Mesh
    EXECUTABLE_ARGS -i RainEvent30.asc -o ${Data_BINARY_DIR}/FileConverter/RainEvent30-materials.vtu -e quad -p materials
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    DIFF_DATA
    RainEvent30-materials.vtu RainEvent30-materials.vtu MaterialIDs MaterialIDs 0 0
)

AddTest(
    NAME Raster2Mesh_Scalars_Test
    PATH FileConverter
    EXECUTABLE Raster2Mesh
    EXECUTABLE_ARGS -i RainEvent30.asc -o ${Data_BINARY_DIR}/FileConverter/RainEvent30-scalars.vtu -e tri -p scalar -n ScalarValues
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    DIFF_DATA
    RainEvent30-scalars.vtu RainEvent30-scalars.vtu ScalarValues ScalarValues 0 0
)
