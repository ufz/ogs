if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1D1bt/m2_1D1bt.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1D2bt/m2_1D2bt.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1Dcreep/m2_1Dcreep.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1Dlozenge/m2_1Dlozenge.prj RUNTIME 5)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1Dlozengebt/m2_1Dlozengebt.prj RUNTIME 29)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1Drelax/m2_1Drelax.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_2Dload/m2_2Dload.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_2Dload/m2_2Dload_ym45.prj RUNTIME 16)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_2Dloadbt/m2_2Dloadbt.prj RUNTIME 25)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_3Dload/m2_3Dload.prj RUNTIME 9)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_3Dloadbt/m2_3Dloadbt.prj RUNTIME 26)
    OgsTest(PROJECTFILE ThermoMechanics/InitialStates/into_initial_state.prj)
    OgsTest(PROJECTFILE ThermoMechanics/InitialStates/equilibrium_restart.prj)
    OgsTest(PROJECTFILE ThermoMechanics/InitialStates/non_equilibrium_initial_state.prj)
    OgsTest(PROJECTFILE ThermoMechanics/EmbeddedAnchorSourceTerm/prism_hex.prj)
    NotebookTest(NOTEBOOKFILE ThermoMechanics/CreepBGRa/CreepAfterExcavation/creep-after-excavation.py RUNTIME 128)
    # Staggered Scheme
    OgsTest(PROJECTFILE ThermoMechanics/StaggeredScheme/TM_Quad/iglu_quarter_plane_strain_quad.prj RUNTIME 11)
    OgsTest(PROJECTFILE ThermoMechanics/StaggeredScheme/CreepAfterExcavation/CreepAfterExcavation.prj RUNTIME 7)
    OgsTest(PROJECTFILE ThermoMechanics/tm1_1Dbeam/tm1_1Dbeam.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoMechanics/tm1_1Dfixa/tm1_1Dfixa.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoMechanics/tm1_1Dfixb/tm1_1Dfixb.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoMechanics/tm1_2Dbeam/tm1_2Dbeam.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoMechanics/tm1_2Dsquare/tm1_2Dsquare.prj RUNTIME 1)
endif()

if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoMechanics/tm1_3Dcube/tm1_3Dcube.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoMechanics/tm1_3Dgravity/tm1_3Dgravity.prj
            RUNTIME 1
    )
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoMechanics/tm1_3Dorigin/tm1_3Dorigin.prj RUNTIME 2)
    OgsTest(PROJECTFILE ThermoMechanics/tm1_3Dsquare/tm1_3Dsquare.prj
            RUNTIME 55
    )
    OgsTest(PROJECTFILE ThermoMechanics/tm2_1D1bt/tm2_1D1bt.prj RUNTIME 10)
    OgsTest(PROJECTFILE ThermoMechanics/tm2_1Dfixc/tm2_1Dfixc.prj RUNTIME 25)
    #--
    OgsTest(PROJECTFILE ThermoMechanics/cube_1e3.prj RUNTIME 2)
    OgsTest(PROJECTFILE ThermoMechanics/iglu_quarter_plane_strain.prj RUNTIME 7)
    OgsTest(PROJECTFILE ThermoMechanics/iglu_axisymmetric_plane_strain.prj)
    OgsTest(PROJECTFILE ThermoMechanics/iglu_quarter_plane_strain_quad.prj RUNTIME 45)
    OgsTest(PROJECTFILE ThermoMechanics/iglu_axisymmetric_plane_strain_quad.prj)
endif()

if(NOT (OGS_USE_LIS OR OGS_USE_MPI OR ENABLE_ASAN))
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/SimpleAxisymmetricCreep/SimpleAxisymmetricCreep.prj RUNTIME 4)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(
        PROJECTFILE ThermoMechanics/CreepBGRa/SimpleAxisymmetricCreep/SimpleAxisymmetricCreepWithAnalyticSolution.prj
        RUNTIME 8
    )
endif()

if(OGS_USE_MFRONT AND NOT (OGS_USE_LIS OR OGS_USE_MPI))
    OgsTest(
        PROJECTFILE ThermoMechanics/CreepBGRa/SimpleAxisymmetricCreep/SimpleAxisymmetricCreepWithAnalyticSolutionMFront.prj
        RUNTIME 6
    )
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/CreepAfterExcavation/CreepAfterExcavation.prj RUNTIME 8)
endif()

# Basic test that MFront models work for TM.
# Linear elastic, no internal state variables, but external temperature.
if(OGS_USE_MFRONT AND NOT (OGS_USE_LIS OR OGS_USE_MPI))
    OgsTest(PROJECTFILE ThermoMechanics/LinearMFront/cube_1e0_lin.prj)
endif()

# Test of a creep law.
if(OGS_USE_MFRONT AND NOT (OGS_USE_LIS OR OGS_USE_MPI))
    OgsTest(PROJECTFILE ThermoMechanics/BDT/cube_1e0_bdt.prj)
endif()
