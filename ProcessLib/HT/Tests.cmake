# Monolithic scheme
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/HT/ConstViscosity/square_5500x5500.prj
            RUNTIME 27
    )
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/IsothermalFluidFlow.prj)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/PressureDiffusionTemperatureDiffusion.prj)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/IsothermalFluidFlowWithGravity.prj)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/PressureParabolicTemperatureParabolic.prj)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/CoupledPressureParabolicTemperatureParabolic.prj)
endif()

# Staggered + XDMF
if(TARGET xdmfdiff AND NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(
        PROJECTFILE
            Parabolic/HT/SimpleSynthetics/XDMF/CoupledPressureParabolicTemperatureParabolicStaggered.prj
    )
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/calculatesurfaceflux_ht_cube_1e3.prj)
    # TODO: relaxed tol. after !5697, maybe regenerate ref. results?
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/calculatesurfaceflux_ht_cube_1e4.prj RUNTIME 73)
    # Staggered scheme
    OgsTest(PROJECTFILE Parabolic/HT/StaggeredCoupling/ADecovalexTHMCBasedHTExample/th_decovalex.prj RUNTIME 2)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/IsothermalFluidFlowStaggered.prj)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/PressureDiffusionTemperatureDiffusionStaggered.prj)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/IsothermalFluidFlowWithGravityStaggered.prj)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/PressureParabolicTemperatureParabolicStaggered.prj)
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/CoupledPressureParabolicTemperatureParabolicStaggered.prj)
    OgsTest(PROJECTFILE Parabolic/HT/FaultedCube/Ra_795_fault_bcgs_jacobi.prj
            RUNTIME 12
    )
    # generateInvalidMediaForHT.py logic moved to PythonSetup.cmake

    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/deactivated_subdomain/HT_DeactivatedSubdomain.prj)
    OgsTest(PROJECTFILE Parabolic/HT/LowerDimensionalFracture/2D_single_fracture_HT.prj RUNTIME 11)
    OgsTest(PROJECTFILE Parabolic/HT/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ThermalDiffusion/TemperatureField.prj RUNTIME 15)
    OgsTest(PROJECTFILE Parabolic/HT/StaggeredCoupling/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow.prj)
    OgsTest(
        PROJECTFILE Parabolic/HT/StaggeredCoupling/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow_via_sub-coupling.xml
        RUNTIME 1
    )
    OgsTest(PROJECTFILE Parabolic/HT/ClassicalTransportExample/classical_transport_example.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/HT/ClassicalTransportExample/classical_transport_example_full_upwind.prj RUNTIME 1)
    OgsTest(
        PROJECTFILE Parabolic/HT/ClassicalTransportExample/classical_transport_example_full_upwind_staggered.prj
        RUNTIME 1
    )
endif()

#MPI/PETSc
if(OGS_USE_MPI)
    OgsTest(
        PROJECTFILE
            Parabolic/HT/StaggeredCoupling/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow_mpi.xml
        WRAPPER mpirun -np 3
        RUNTIME 2
    )
endif()

if(OGS_USE_PIP AND NOT (OGS_USE_MPI OR OGS_USE_LIS))
    if(NOT EXISTS ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/HT_specific_heat_capacity_viscosity_porosity.prj)
        execute_process(
            COMMAND
                uv run python
                ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/generateInvalidMediaForHT.py
            WORKING_DIRECTORY
                ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles
            RESULT_VARIABLE GEN_INVALID_RES
        )
        if(GEN_INVALID_RES EQUAL 0)
            message(STATUS "generateInvalidMediaForHT.py succeeded.")
        else()
            message(SEND_ERROR "generateInvalidMediaForHT.py failed with status ${GEN_INVALID_RES}.")
        endif()
    endif()
    file(GLOB HT_INVALID_PRJ_FILES
            ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/*.prj
    )
    foreach(ht_invalid_prj_file ${HT_INVALID_PRJ_FILES})
        file(
            RELATIVE_PATH ht_invalid_prj_file_rel
            ${Data_SOURCE_DIR}
            ${ht_invalid_prj_file}
        )
        OgsTest(
            PROJECTFILE ${ht_invalid_prj_file_rel}
            RUNTIME 1
            PROPERTIES WILL_FAIL TRUE
        )
    endforeach()
endif()
