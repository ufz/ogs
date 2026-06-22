# Liquid flow
if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/LineDirichletNeumannBC/line_dirichlet_neumannBC.prj)
    OgsTest(
        PROJECTFILE Parabolic/LiquidFlow/PressureBCatCornerOfAnisotropicSquare/pressureBC_at_corner_of_anisotropic_square.prj
    )
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/DrainageExcavation/drainage_LiquidFlow.prj)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/GravityDriven/gravity_driven.prj)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/GravityDriven/gravity_driven_XZ.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/AxiSymTheis/axisym_theis.prj)
endif()

if(OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/AxiSymTheis/axisym_theis.prj WRAPPER mpirun -np 1)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/BuildupTest/buildup_test.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/GravityDriven3D/anisotropic_gravity_driven3D.prj RUNTIME 53)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/GravityDriven3D/isotropic_gravity_driven3D.prj RUNTIME 51)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/TimeIntervalDirichletBC/TimeIntervalDirichletBC.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Verification/h1_1Dsource/h1_1Dsource.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Verification/h1_1Dsteady/h1_1Dsteady.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Verification/h1_3Dhydstat/h1_3Dhydstat.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Verification/h2_1D1bt/h2_1D1bt.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Verification/h2_1D2bt/h2_1D2bt.prj)
endif()

#===============================================================================
# PETSc/MPI
if(OGS_USE_MPI)
    OgsTest(
        PROJECTFILE Parabolic/LiquidFlow/LineDirichletNeumannBC/line_dirichlet_neumannBC_petsc.xml
        WRAPPER mpirun -np 1
    )
    OgsTest(
        PROJECTFILE Parabolic/LiquidFlow/GravityDriven/gravity_driven_petsc.xml
        WRAPPER mpirun -np 1
    )
    OgsTest(
        PROJECTFILE Parabolic/LiquidFlow/PressureBCatCornerOfAnisotropicSquare/pressureBC_at_corner_of_anisotropic_square.prj
        WRAPPER mpirun -np 1
    )
    OgsTest(
        PROJECTFILE Parabolic/LiquidFlow/GravityDriven3D/anisotropic_gravity_driven3D.prj
        WRAPPER mpirun -np 1
        RUNTIME 115
    )
    OgsTest(
        PROJECTFILE Parabolic/LiquidFlow/GravityDriven3D/isotropic_gravity_driven3D.prj
        WRAPPER mpirun -np 1
        RUNTIME 130
    )
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/TimeIntervalDirichletBC/TimeIntervalDirichletBC.prj WRAPPER mpirun -np 1)
endif()

# Dupuit
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Unconfined_Aquifer/BC_BC/TestSet_01.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Unconfined_Aquifer/BC_BC_RECHARGE/TestSet_01.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Unconfined_Aquifer/BC_BC_RECHARGE2/TestSet_01.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Unconfined_Aquifer/BC_BC_STORAGE/TestSet_01.prj)
    OgsTest(
        PROJECTFILE
            Parabolic/LiquidFlow/TimeDependentHeterogeneousBoundaryConditions/TimeDependentHeterogeneousBoundaryConditions.prj
    )
    OgsTest(
        PROJECTFILE
            Parabolic/LiquidFlow/TimeDependentHeterogeneousSourceTerm/TimeDependentHeterogeneousSourceTerm.prj
    )
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Flux/3D/Pyramid/cuboid_1x1x1_pyramid_6000_calculatesurfaceflux.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Flux/2D/square_1e1_calculatesurfaceflux.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Flux/2D/square_1.8e1_calculatesurfaceflux.prj)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Flux/3D/Hex/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_multiple_outputs.xml)

    AddTest(
        NAME LiquidFlow_Flux_3D_HEX_MultipleOutputs_xdmf
        PATH Parabolic/LiquidFlow/Flux/3D/Hex
        EXECUTABLE ogs
        EXECUTABLE_ARGS
            cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_multiple_outputs.xml
        WRAPPER time
        REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
        TESTER xdmfdiff
        DIFF_DATA
            top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_cuboid_1x1x1_hex_27.xdmf
            top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_cuboid_1x1x1_hex_27.xdmf
            pressure pressure 1e-10 1e-15 0 0
            top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_cuboid_1x1x1_hex_27_bottom_boundary.xdmf
            top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_cuboid_1x1x1_hex_27_bottom_boundary.xdmf
            pressure pressure 1e-7 1e-13 0 0
            top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_cuboid_1x1x1_hex_27_top_boundary.xdmf
            top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_cuboid_1x1x1_hex_27_top_boundary.xdmf
            pressure pressure 1e-7 1e-13 0 0
        PROPERTIES DEPENDS ogs-Parabolic/LiquidFlow/Flux/3D/Hex/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_multiple_outputs
    )

    if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
        OgsTest(
            PROJECTFILE Parabolic/LiquidFlow/Flux/3D/Hex/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_invalid_multiple_outputs.xml
            NO_TEST_DEFINITION
            PROPERTIES PASS_REGULAR_EXPRESSION "Output configuration paths are not unique. This will lead to overwritten results or invalid"
        )
    endif()

    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Flux/3D/Hex/MultipleOutputsDifferentVariablesSections/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_multiple_outputs_different_variables.xml)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/Flux/3D/Hex/MultipleOutputsDifferentVariablesSections/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_empty_output_timesteps_only_fixed_output.xml)
    set(EXPECTED_FILES
        empty_output_timesteps_only_fixed_output_times_config_cuboid_1x1x1_hex_27_ts_0_t_0.000000.vtu
        empty_output_timesteps_only_fixed_output_times_config_cuboid_1x1x1_hex_27_ts_5_t_43200.000000.vtu
        empty_output_timesteps_only_fixed_output_times_config_cuboid_1x1x1_hex_27_ts_10_t_86400.000000.vtu
    )
    set(DIR_TO_CHECK
        ${PROJECT_BINARY_DIR}/Tests/Data/Parabolic/LiquidFlow/Flux/3D/Hex/MultipleOutputsDifferentVariablesSections_52a2ff3d
    )
    set(FILE_PREFIX
        empty_output_timesteps_only_fixed_output_times_config_cuboid_1x1x1_hex_27_ts_
    )
    add_test(
        NAME check_files-ogs-Parabolic/LiquidFlow/Flux/3D/Hex
        COMMAND
            ${CMAKE_COMMAND} "-DEXPECTED_FILES=${EXPECTED_FILES}"
            -DFILE_PREFIX=${FILE_PREFIX} -DDIR_TO_CHECK=${DIR_TO_CHECK} -P
            ${PROJECT_SOURCE_DIR}/scripts/cmake/test/CheckCreatedFiles.cmake
    )
    set_tests_properties(
        check_files-ogs-Parabolic/LiquidFlow/Flux/3D/Hex
        PROPERTIES
            DEPENDS
            ogs-Parabolic/LiquidFlow/Flux/3D/Hex/MultipleOutputsDifferentVariablesSections/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_empty_output_timesteps_only_fixed_output
    )
endif()

if(OGS_USE_MPI)
    OgsTest(
        PROJECTFILE
            Parabolic/LiquidFlow/Flux/3D/Hex/Parallel/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet.prj
        WRAPPER mpirun -np 2
    )
endif()

AddTest(
    NAME LiquidFlow_SimpleSynthetics_XDMF
    PATH Parabolic/LiquidFlow/SimpleSynthetics/XDMF
    EXECUTABLE ogs
    EXECUTABLE_ARGS FunctionParameterTest_XDMF.prj
    WRAPPER time
    TESTER xdmfdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    DIFF_DATA
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf pressure pressure 1e-7 1e-13 1 1
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf VolumetricFlowRate VolumetricFlowRate 1e-7 1e-13 2 2
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf MaterialIDs MaterialIDs 1e-7 1e-13 3 3
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf v v 1e-7 1e-13 4 4
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf pressure pressure 1e-7 1e-13 2 2
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf bulk_element_ids bulk_element_ids 1e-7 1e-13 3 3
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf bulk_node_ids bulk_node_ids 1e-7 1e-13 2 2
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf pressure pressure 1e-7 1e-13 4 4
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf bulk_element_ids bulk_element_ids 1e-7 1e-13 3 3
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf bulk_node_ids bulk_node_ids 1e-7 1e-13 2 2
)

AddTest(
    NAME LiquidFlow_SimpleSynthetics_XDMF_compression_off
    PATH Parabolic/LiquidFlow/SimpleSynthetics/XDMF_compression_off
    EXECUTABLE ogs
    EXECUTABLE_ARGS FunctionParameterTest_XDMF.prj
    WRAPPER time
    TESTER xdmfdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    DIFF_DATA
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf pressure pressure 1e-7 1e-13 1 1
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf VolumetricFlowRate VolumetricFlowRate 1e-7 1e-13 2 2
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf MaterialIDs MaterialIDs 1e-7 1e-13 3 3
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf v v 1e-7 1e-13 1 1
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf pressure pressure 1e-7 1e-13 2 2
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf bulk_element_ids bulk_element_ids 1e-7 1e-13 0 0
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf bulk_node_ids bulk_node_ids 1e-7 1e-13 1 1
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf pressure pressure 1e-7 1e-13 2 2
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf bulk_element_ids bulk_element_ids 1e-7 1e-13 0 0
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf bulk_node_ids bulk_node_ids 1e-7 1e-13 1 1
)

# Single mesh XDMF MPI tests
foreach(mesh bulk left right top bottom)
    AddTest(
        NAME LiquidFlow_SimpleSynthetics_XDMF_MPI_${mesh}
        PATH Parabolic/LiquidFlow/SimpleSynthetics/XDMF_MPI/3/${mesh}
        EXECUTABLE ogs
        EXECUTABLE_ARGS ${mesh}.xml -m ../
        WRAPPER mpirun
        WRAPPER_ARGS -np 3
        TESTER xdmfdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        ${mesh}_${mesh}.xdmf  ${mesh}_${mesh}.xdmf  pressure pressure 1e-12 1e-12 1 1
        ${mesh}_${mesh}.xdmf  ${mesh}_${mesh}.xdmf  v v 1e-12 1e-12 1 1
        ${mesh}_${mesh}.xdmf  ${mesh}_${mesh}.xdmf  VolumetricFlowRate VolumetricFlowRate 1e-12 1e-12 1 1
    )
endforeach()

# Two mesh XDMF MPI tests
foreach(mesh1 bulk left right top bottom)
    foreach(mesh2 bulk left right top bottom)
        if(NOT "${mesh1}" STREQUAL "${mesh2}")
            AddTest(
                NAME LiquidFlow_SimpleSynthetics_XDMF_MPI_${mesh1}_${mesh2}
                PATH Parabolic/LiquidFlow/SimpleSynthetics/XDMF_MPI/3/${mesh1}_${mesh2}
                EXECUTABLE ogs
                EXECUTABLE_ARGS ${mesh1}_${mesh2}.xml -m ../
                WRAPPER mpirun
                WRAPPER_ARGS -np 3
                TESTER xdmfdiff
                REQUIREMENTS OGS_USE_MPI
                DIFF_DATA
                ${mesh1}_${mesh1}.xdmf  ${mesh1}_${mesh1}.xdmf  pressure pressure 1e-12 1e-12 1 1
                ${mesh1}_${mesh1}.xdmf  ${mesh1}_${mesh1}.xdmf  v v 1e-12 1e-12 1 1
                ${mesh1}_${mesh1}.xdmf  ${mesh1}_${mesh1}.xdmf  VolumetricFlowRate VolumetricFlowRate 1e-12 1e-12 1 1
                ${mesh1}_${mesh2}.xdmf  ${mesh1}_${mesh2}.xdmf  pressure pressure 1e-12 1e-12 1 1
                ${mesh1}_${mesh2}.xdmf  ${mesh1}_${mesh2}.xdmf  v v 1e-12 1e-12 1 1
                ${mesh1}_${mesh2}.xdmf  ${mesh1}_${mesh2}.xdmf  VolumetricFlowRate VolumetricFlowRate 1e-12 1e-12 1 1
            )
        endif()
    endforeach()
endforeach()

# Three mesh XDMF MPI tests
foreach(mesh1 bulk left right top bottom)
    foreach(mesh2 bulk left right top bottom)
        foreach(mesh3 bulk left right top bottom)
            if(NOT "${mesh1}" STREQUAL "${mesh2}" AND NOT "${mesh1}" STREQUAL "${mesh3}" AND NOT "${mesh2}" STREQUAL "${mesh3}")
                AddTest(
                    NAME LiquidFlow_SimpleSynthetics_XDMF_MPI_${mesh1}_${mesh2}_${mesh3}
                    PATH Parabolic/LiquidFlow/SimpleSynthetics/XDMF_MPI/3/${mesh1}_${mesh2}_${mesh3}
                    EXECUTABLE ogs
                    EXECUTABLE_ARGS ${mesh1}_${mesh2}_${mesh3}.xml -m ../
                    WRAPPER mpirun
                    WRAPPER_ARGS -np 3
                    TESTER xdmfdiff
                    REQUIREMENTS OGS_USE_MPI
                    DIFF_DATA
                    ${mesh1}_${mesh1}.xdmf  ${mesh1}_${mesh1}.xdmf  pressure pressure 1e-12 1e-12 1 1
                    ${mesh1}_${mesh1}.xdmf  ${mesh1}_${mesh1}.xdmf  v v 1e-12 1e-12 1 1
                    ${mesh1}_${mesh1}.xdmf  ${mesh1}_${mesh1}.xdmf  VolumetricFlowRate VolumetricFlowRate 1e-12 1e-12 1 1
                    ${mesh1}_${mesh2}.xdmf  ${mesh1}_${mesh2}.xdmf  pressure pressure 1e-12 1e-12 1 1
                    ${mesh1}_${mesh2}.xdmf  ${mesh1}_${mesh2}.xdmf  v v 1e-12 1e-12 1 1
                    ${mesh1}_${mesh2}.xdmf  ${mesh1}_${mesh2}.xdmf  VolumetricFlowRate VolumetricFlowRate 1e-12 1e-12 1 1
                    ${mesh1}_${mesh3}.xdmf  ${mesh1}_${mesh3}.xdmf  pressure pressure 1e-12 1e-12 1 1
                    ${mesh1}_${mesh3}.xdmf  ${mesh1}_${mesh3}.xdmf  v v 1e-12 1e-12 1 1
                    ${mesh1}_${mesh3}.xdmf  ${mesh1}_${mesh3}.xdmf  VolumetricFlowRate VolumetricFlowRate 1e-12 1e-12 1 1
                )
            endif()
        endforeach()
    endforeach()
endforeach()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/PrimaryVariableConstraintDirichletBC/cuboid_1x1x1_hex_1000_Dirichlet_Dirichlet_1.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/PrimaryVariableConstraintDirichletBC/cuboid_1x1x1_hex_1000_Dirichlet_Dirichlet_2.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/PrimaryVariableConstraintDirichletBC/cuboid_1x1x1_hex_1000_Dirichlet_Dirichlet_3.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/FunctionParameterTest.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/BlockingConductingFracture/block_conduct_frac.prj)
    NotebookTest(NOTEBOOKFILE Parabolic/LiquidFlow/BlockingConductingFracture/BlockingConductingFracture.py RUNTIME 9)
    NotebookTest(NOTEBOOKFILE Parabolic/LiquidFlow/roughFracture/roughFracture_benchmark.py RUNTIME 600)
    # inclined mesh
    OgsTest(
        PROJECTFILE
            Parabolic/LiquidFlow/InclinedMeshElements/Inclined2DMesh/hydrostatic_flow_in_inclined_2D_plane.prj
    )
    OgsTest(
        PROJECTFILE
            Parabolic/LiquidFlow/InclinedMeshElements/Inclined2DMesh/transient_flow_in_inclined_2D_plane.prj
    )
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/InclinedMeshElements/FractureIn3D/fractures_in_3D.prj RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/InclinedMeshElements/1Din3D/line_fractures_in_3D.prj RUNTIME 1)
    OgsTest(PROJECTFILE Utils/GMSH2OGS/quadratic_mesh_assembly_test.prj
            RUNTIME 1
    )
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/RasterParameter/CoarseRasterHomogeneous.xml RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/RasterParameter/FineRasterHomogeneous.xml RUNTIME 1)
    # fine raster and coarse raster input should produce exactly the same output
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/RasterParameter/CoarseRasterHeterogeneous.xml RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/RasterParameter/FineRasterHeterogeneous.xml RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/GasFlow/gas_flow.prj RUNTIME 1)
endif()
