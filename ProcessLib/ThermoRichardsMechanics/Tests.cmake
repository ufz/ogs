if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoRichardsMechanics/LinearMechanics/mechanics_linear.prj)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/FullySaturatedFlowMechanics/flow_fully_saturated.prj)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/RichardsFlow2D/RichardsFlow_2d_small.prj RUNTIME 3)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/A2/A2.prj RUNTIME 8)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/Porosity/deformation_dependent_porosity.prj RUNTIME 3)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/Porosity/deformation_dependent_porosity_swelling.prj RUNTIME 6)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/Porosity/deformation_temperature_dependent_porosity_swelling.prj RUNTIME 11)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/OrthotropicSwelling/orthotropic_swelling_xy.xml)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/OrthotropicSwelling/orthotropic_swelling_xyz.xml)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/anisotropic_thermal_expansion/aniso_expansion.prj RUNTIME 1)
    OgsTest(
        PROJECTFILE ThermoRichardsMechanics/anisotropic_thermal_expansion/aniso_expansion_expansivity_matrix.xml
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE ThermoRichardsMechanics/anisotropic_thermal_expansion/aniso_expansion_expansivity_matrix_z90.xml
        RUNTIME 1
    )
    OgsTest(PROJECTFILE ThermoRichardsMechanics/anisotropic_thermal_expansion/aniso_expansion_x45.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/anisotropic_thermal_expansion/aniso_expansion_y45.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/anisotropic_thermal_expansion/aniso_expansion_z45.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/LiakopoulosHM/liakopoulos.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/LiakopoulosHM/liakopoulos_restart.xml RUNTIME 1)
endif()

if(OGS_USE_PETSC)
    OgsTest(
        PROJECTFILE ThermoRichardsMechanics/LiakopoulosPETSc/liakopoulos_mixElem_mumps.prj
        WRAPPER mpirun -np 2
        RUNTIME 2
        LABELS "petsc-mumps"
        NAME_SUFFIX LiakopoulosMixedElementsPETSc
    )
endif()

if(NOT OGS_USE_MPI)
    OgsTest(
        PROJECTFILE
            ThermoRichardsMechanics/Simple3DThermoMechanicsFromTM/cube_1e3_stress_analysis.xml
        RUNTIME 5
    )
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(
        PROJECTFILE ThermoRichardsMechanics/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow.prj
        RUNTIME 1
    )
    OgsTest(PROJECTFILE ThermoRichardsMechanics/PointHeatSource/point_heat_source_2D.prj RUNTIME 4)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/TaskCDECOVALEX2023/Decovalex-0.prj RUNTIME 15)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/CTF1/CTF1.prj RUNTIME 3)
endif()
if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/SimpleAxisymmetricCreep/SimpleAxisymmetricCreep.prj RUNTIME 6)
endif()

#PETSc
if(OGS_USE_MPI)
    OgsTest(
        PROJECTFILE
            ThermoRichardsMechanics/Simple3DThermoMechanicsFromTM/cube_1e3.prj
        WRAPPER mpirun -np 3
        RUNTIME 40 LABELS "petsc-mumps"
    )
    OgsTest(
        PROJECTFILE
            ThermoRichardsMechanics/FullySaturatedFlowMechanics/PETSc/flow_fully_saturated_petsc.prj
        WRAPPER mpirun -np 2
        RUNTIME 10
    )
    OgsTest(
        PROJECTFILE ThermoRichardsMechanics/PointHeatSource/point_heat_source_2D.prj
        PATCH_FILES point_heat_source_2D_non_submesh_r_output.xml
        NAME_SUFFIX non_submesh_r_output
        WRAPPER mpirun -np 3
        RUNTIME 5
        LABELS "petsc-mumps"
    )
    OgsTest(
        PROJECTFILE
            ThermoRichardsMechanics/PointHeatSource/point_heat_source_2D_gml.prj
        WRAPPER mpirun -np 3
        RUNTIME 5 LABELS "petsc-mumps"
    )
    OgsTest(
        PROJECTFILE
            ThermoRichardsMechanics/TaskCDECOVALEX2023/Decovalex-0_mpi.xml
        WRAPPER mpirun -np 3
        RUNTIME 10 LABELS "petsc-mumps"
    )
endif()
# ThermoRichardsMechanics; thermo_osmosis and thermo_filtration effects, linear poroelastic, column consolidation
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoRichardsMechanics/ThermoOsmosis/Column.prj RUNTIME 9)
    # ThermoRichardsMechanics; test for removing body force from displacement equation
    OgsTest(PROJECTFILE ThermoRichardsMechanics/BodyForce/square.prj RUNTIME 1)
endif()

if(NOT OGS_USE_LIS)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/BodyForce/square_total_stress_test.xml RUNTIME 1)
    OgsTest(
        PROJECTFILE ThermoRichardsMechanics/Simple3DThermoMechanicsFromTM/cube_1e3.prj
        PATCH_FILES 3D_axially_symmetric_fail_test_patch.xml
        NAME_SUFFIX axially_symmetric_fail
        PROPERTIES
            PASS_REGULAR_EXPRESSION "3D mesh cannot be axially symmetric."
    )

endif()

if(OGS_USE_MFRONT AND (NOT OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MultiMaterialEhlers/square_1e1_2_matIDs.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MultiMaterialEhlers/square_1e1_2_matIDs_restart.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/ThermoPoroElasticity/uniaxial_isothermal_drainage_imbibition_basic_mfront_model_ctest.xml)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/ThermoPoroElasticity/uniaxial_isothermal_drainage_imbibition_extended_mfront_model_ctest.xml)

    if (NOT OGS_USE_MPI)
        OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/A2/A2.xml RUNTIME 8)
        if(NOT OGS_COVERAGE)
            OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/A2/A2_effective_stress0.xml RUNTIME 8)
        endif()
    endif()

    add_subdirectory(
        "${CMAKE_SOURCE_DIR}/Tests/Data/ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/MFrontBehaviour/"
        "${CMAKE_BINARY_DIR}/Tests/Data/ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/MFrontBehaviour/"
    )

    OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/0d_confined_compression/confined_compression.prj RUNTIME 19)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/0d_resaturation/resaturation.prj RUNTIME 4)

    if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
        # Disabled for PETSc, because large numerical errors are observed occasionally.
        # The model is tested sufficiently in other configurations, and other tests are run with PETSc, too.
        OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/1d_column_resaturation/bentonite_column.prj RUNTIME 1500)
    endif()

    OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/1d_column_restart/bentonite_column_restart.xml RUNTIME 5)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/1d_column_restart/bentonite_column_restart_fail.xml PROPERTIES PASS_REGULAR_EXPRESSION "Absolute and relative error [(]maximum norm[)] are larger than the corresponding thresholds 8[.]000000000000000e-03 and 2[.]000000000000000e-02[.]
" RUNTIME 10)
endif()

if (NOT (OGS_USE_MPI OR OGS_USE_LIS))
    if(NOT WIN32) # TODO: Remove after 6.5.5 release and fix benchmark on win
        OgsTest(PROJECTFILE ThermoRichardsMechanics/Mockup2D/mockup.prj RUNTIME 15)
    endif()
    OgsTest(PROJECTFILE ThermoRichardsMechanics/Mockup2D/mockup_restart.xml RUNTIME 5)
endif()
