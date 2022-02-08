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

foreach(process ${_processes_list})
    option(OGS_BUILD_PROCESS_${process} "Build process ${process}" ON)
    if(OGS_BUILD_PROCESS_${process})
        list(APPEND _enabled_processes ${process})
    endif()
endforeach()

if(NOT OGS_BUILD_CLI)
    set(OGS_BUILD_PROCESSES OFF)
    message(
        STATUS
            "ATTENTION: OGS_BUILD_CLI=OFF -> OGS_BUILD_PROCESSES is set to OFF too.\n"
            "   If cli is switched on again, remember to switch processes back to on \n"
            "   too with -DOGS_BUILD_PROCESSES=ON!"
    )
endif()

if(DEFINED CACHE{OGS_BUILD_PROCESSES})
    if("${OGS_BUILD_PROCESSES}" STREQUAL "ON" OR "${OGS_BUILD_PROCESSES}"
                                                 STREQUAL "OFF"
    )
        if(OGS_BUILD_PROCESSES)
            set(_enabled_processes ${_processes_list})
        else()
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
    endif()
    foreach(process ${_processes_list})
        if(${process} IN_LIST _enabled_processes)
            set(OGS_BUILD_PROCESS_${process} ON CACHE BOOL "" FORCE)
        else()
            set(OGS_BUILD_PROCESS_${process} OFF CACHE BOOL "" FORCE)
        endif()
    endforeach()
    unset(OGS_BUILD_PROCESSES CACHE)
endif()

if("RichardsComponentTransport" IN_LIST _enabled_processes
   OR "ComponentTransport" IN_LIST _enabled_processes
)
    set(_build_chemistry_lib ON)
endif()

# Print summary of enabled processes.
function(printEnabledProcesses)
    list(LENGTH _enabled_processes num_enabled_processes)
    list(LENGTH _processes_list num_processes_list)
    if(${num_enabled_processes} EQUAL 0)
        message(STATUS "All processes have been disabled!\n")
        return()
    endif()
    if(${num_processes_list} EQUAL ${num_enabled_processes})
        message(STATUS "All processes have been enabled.\n")
        return()
    endif()
    message(STATUS "The following processes have been enabled:\n")
    foreach(process ${_enabled_processes})
        message(" * ${process}")
    endforeach()
    message("")
endfunction()
