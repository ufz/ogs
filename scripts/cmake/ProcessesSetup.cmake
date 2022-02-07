# A list of processes to be build. Also used in the ProcessLib to select
# processes to be build.
set(_processes_list
    ComponentTransport
    StokesFlow
    HT
    HeatConduction
    HeatTransportBHE
    HydroMechanics
    LiquidFlow
    LIE
    ThermoRichardsMechanics
    PhaseField
    RichardsComponentTransport
    RichardsFlow
    RichardsMechanics
    SmallDeformation
    SmallDeformationNonlocal
    SteadyStateDiffusion
    TES
    TH2M
    ThermalTwoPhaseFlowWithPP
    ThermoHydroMechanics
    ThermoMechanicalPhaseField
    ThermoMechanics
    ThermoRichardsFlow
    TwoPhaseFlowWithPP
    TwoPhaseFlowWithPrho
)

set(OGS_BUILD_PROCESSES ON
    CACHE STRING "Semicolon-separated list of processes to build"
)

# Drop-down (cmake-gui), cycle-through list (ccmake) for easiere interactive
# selection. For enabling multiple selected processes cmake usage via cli is
# required!
set_property(
    CACHE OGS_BUILD_PROCESSES PROPERTY STRINGS ON OFF ${_processes_list}
)

if(NOT OGS_BUILD_CLI)
    set(OGS_BUILD_PROCESSES OFF "" CACHE STRING "" FORCE)
    message(
        STATUS
            "ATTENTION: OGS_BUILD_CLI=OFF -> OGS_BUILD_PROCESSES is set to OFF too.\n"
            "   If cli is switched on again, remember to switch processes back to on \n"
            "   too with -DOGS_BUILD_PROCESSES=ON!"
    )
endif()

if("${OGS_BUILD_PROCESSES}" STREQUAL "ON" OR "${OGS_BUILD_PROCESSES}" STREQUAL
                                             "OFF"
)
    if(OGS_BUILD_PROCESSES)
        message(STATUS "All processes enabled.")
        set(_enabled_processes ${_processes_list})
    else()
        message(STATUS "All processes disabled.")
        unset(_enabled_processes)
    endif()
else()
    foreach(process ${OGS_BUILD_PROCESSES})
        if(NOT "${process}" IN_LIST _processes_list)
            message(
                FATAL_ERROR
                    "${process} given in OGS_BUILD_PROCESSES is "
                    "not a valid process name! Valid names are ${_processes_list}"
            )
        endif()
    endforeach()
    set(_enabled_processes ${OGS_BUILD_PROCESSES})
    message(STATUS "Enabled processes: ${OGS_BUILD_PROCESSES}")
endif()
