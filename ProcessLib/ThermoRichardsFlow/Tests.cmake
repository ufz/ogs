if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoRichardsFlow/HT/SimpleSynthetics/PressureDiffusionTemperatureDiffusion.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsFlow/HT/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow.prj RUNTIME 2)
    OgsTest(PROJECTFILE ThermoRichardsFlow/RichardsFlow2D/RichardsFlow_2d_small.prj RUNTIME 9)
    OgsTest(PROJECTFILE ThermoRichardsFlow/RichardsFlow2D/RichardsFlow_2d_small_Picard.prj RUNTIME 6)
    OgsTest(PROJECTFILE ThermoRichardsFlow/RichardsFlow2D/RichardsFlow_2d_compare_ogs5.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsFlow/SimplifiedMechanics/TRuni_saturated.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsFlow/SimplifiedMechanics/TRuni_unsaturated.prj RUNTIME 1)
endif()

AddTest(
    NAME ThermoRichardsFlow_HeatTransportInStationaryFlow_compare
    PATH ThermoRichardsFlow/HT/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS <SOURCE_PATH>/HeatTransportInStationaryFlow_compare.xml -l debug
    # auxiliary log files will be put there
    WORKING_DIRECTORY <BUILD_PATH>
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    LABELS NO_PARALLEL_ASSEMBLY
    RUNTIME 17
    DIFF_DATA
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu temperature  temperature 5e-3 1e-8
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu pressure  pressure 5e-3 1e-8
)
# like above but with strict Jacobian comparison tolerances to make the simulation fail
AddTest(
    NAME ThermoRichardsFlow_HeatTransportInStationaryFlow_compare_strict_fail
    PATH ThermoRichardsFlow/HT/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS <SOURCE_PATH>/HeatTransportInStationaryFlow_compare_strict_fail.xml
    # auxiliary log files will be put there
    WORKING_DIRECTORY <BUILD_PATH>
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    LABELS NO_PARALLEL_ASSEMBLY
    PROPERTIES
        PASS_REGULAR_EXPRESSION
        "OGS failed, because the two Jacobian implementations returned different results."
)
# strict Jacobian comparison tolerances but keep going
AddTest(
    NAME ThermoRichardsFlow_HeatTransportInStationaryFlow_compare_strict_nofail
    PATH ThermoRichardsFlow/HT/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS <SOURCE_PATH>/HeatTransportInStationaryFlow_compare_strict_nofail.xml
    # auxiliary log files will be put there
    WORKING_DIRECTORY <BUILD_PATH>
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    LABELS NO_PARALLEL_ASSEMBLY
    RUNTIME 17
    DIFF_DATA
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu temperature  temperature 5e-3 1e-8
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu pressure  pressure 5e-3 1e-8
)

if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoRichardsFlow/SimplifiedMechanics/TRcustom_unsaturated.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsFlow/SimplifiedMechanics/TRhyd_saturated.prj RUNTIME 1)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoRichardsFlow/SimplifiedMechanics/TRhyd_unsaturated.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsFlow/SimplifiedMechanics/TRuni_unsaturated_bishopstest.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsFlow/SimplifiedMechanics/TRhyd_unsaturated_bishopstest.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsFlow/TaskCDECOVALEX2023/Decovalex-0-TRF.prj RUNTIME 7)
    OgsTest(PROJECTFILE ThermoRichardsFlow/TaskCDECOVALEX2023/Decovalex-0-TRF_DeVries.xml RUNTIME 7)
endif()

AddTest(
    NAME ThermoRichardsFlow_TaskCDECOVALEX2023_Picard
    PATH ThermoRichardsFlow/TaskCDECOVALEX2023/WithPicardNonLinearSolverAndPETSc
    EXECUTABLE ogs
    EXECUTABLE_ARGS Decovalex-0-TRF.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 9
    DIFF_DATA
    expected_Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-THuni-0_ts_10_t_864000.000000.vtu pressure pressure 3e5 15
    expected_Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-THuni-0_ts_10_t_864000.000000.vtu saturation saturation 1e-2 2e-3
    expected_Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-THuni-0_ts_10_t_864000.000000.vtu temperature temperature 1e-2 1.0
)

#PETSc
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoRichardsFlow/SimplifiedMechanics/TRuni_unsaturated.prj WRAPPER mpirun -np 1 RUNTIME 5)
    OgsTest(
        PROJECTFILE
            ThermoRichardsFlow/TaskCDECOVALEX2023/WithPicardNonLinearSolverAndPETSc/Decovalex-0-TRF.prj
        WRAPPER mpirun -np 3
        RUNTIME 5
    )
    OgsTest(PROJECTFILE ThermoRichardsFlow/TaskCDECOVALEX2023/Decovalex-0-TRF.prj WRAPPER mpirun -np 1)
endif()
# ThermoRichardsFlow; thermo_osmosis and thermo_filtration effects, linear poroelastic, column consolidation
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoRichardsFlow/ThermoOsmosis/Column.prj RUNTIME 6)
endif()
