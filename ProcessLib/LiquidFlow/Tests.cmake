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

if(TARGET xdmfdiff AND NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(
        PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/XDMF/FunctionParameterTest_XDMF.prj
    )

    OgsTest(
        PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/XDMF/separate_static_dynamic_data_h5_writing.xml
    )

    OgsTest(
        PROJECTFILE
            Parabolic/LiquidFlow/SimpleSynthetics/XDMF_compression_off/FunctionParameterTest_XDMF.prj
    )
endif()

if(OGS_USE_MPI AND TARGET xdmfdiff)
    set(XDMF_MPI_MESH_DIR
        ${Data_SOURCE_DIR}/Parabolic/LiquidFlow/SimpleSynthetics/XDMF_MPI/3
    )

    # Single mesh XDMF MPI tests
    foreach(mesh bulk left right top bottom)
        OgsTest(
            PROJECTFILE
                Parabolic/LiquidFlow/SimpleSynthetics/XDMF_MPI/3/${mesh}/LiquidFlow.prj
            WRAPPER mpirun -np 3
            PATCH_FILES ${mesh}.xml ../xdmfdiff_test_definition.xml
            NAME_SUFFIX ${mesh}
            EXECUTABLE_ARGS -m ${XDMF_MPI_MESH_DIR}
        )
    endforeach()

    # Two mesh XDMF MPI tests
    foreach(mesh1 bulk left right top bottom)
        foreach(mesh2 bulk left right top bottom)
            if(NOT "${mesh1}" STREQUAL "${mesh2}")
                OgsTest(
                    PROJECTFILE
                        Parabolic/LiquidFlow/SimpleSynthetics/XDMF_MPI/3/${mesh1}_${mesh2}/LiquidFlow.prj
                    WRAPPER mpirun -np 3
                    PATCH_FILES ${mesh1}_${mesh2}.xml
                                ../xdmfdiff_test_definition.xml
                    NAME_SUFFIX ${mesh1}_${mesh2}
                    EXECUTABLE_ARGS -m ${XDMF_MPI_MESH_DIR}
                )
            endif()
        endforeach()
    endforeach()

    # Three mesh XDMF MPI tests
    foreach(mesh1 bulk left right top bottom)
        foreach(mesh2 bulk left right top bottom)
            foreach(mesh3 bulk left right top bottom)
                if(NOT "${mesh1}" STREQUAL "${mesh2}" AND NOT "${mesh1}" STREQUAL "${mesh3}" AND NOT "${mesh2}" STREQUAL "${mesh3}")
                    OgsTest(
                        PROJECTFILE
                            Parabolic/LiquidFlow/SimpleSynthetics/XDMF_MPI/3/${mesh1}_${mesh2}_${mesh3}/LiquidFlow.prj
                        WRAPPER mpirun -np 3
                        PATCH_FILES ${mesh1}_${mesh2}_${mesh3}.xml
                                    ../xdmfdiff_test_definition.xml
                        NAME_SUFFIX ${mesh1}_${mesh2}_${mesh3}
                        EXECUTABLE_ARGS -m ${XDMF_MPI_MESH_DIR}
                    )
                endif()
            endforeach()
        endforeach()
    endforeach()
endif()

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
