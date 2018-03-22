include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

add_subdirectory(ComponentTransport)
add_subdirectory(GroundwaterFlow)
add_subdirectory(HT)
add_subdirectory(HeatConduction)
add_subdirectory(HydroMechanics)
add_subdirectory(LiquidFlow)
add_subdirectory(LIE)
add_subdirectory(PhaseField)
add_subdirectory(RichardsComponentTransport)
add_subdirectory(RichardsFlow)
add_subdirectory(SmallDeformation)
add_subdirectory(TES)
add_subdirectory(ThermalTwoPhaseFlowWithPP)
add_subdirectory(ThermoMechanics)
add_subdirectory(TwoPhaseFlowWithPP)
add_subdirectory(TwoPhaseFlowWithPrho)

add_library(Processes INTERFACE)
target_link_libraries(Processes INTERFACE
    ComponentTransport
    GroundwaterFlow
    HT
    HeatConduction
    HydroMechanics
    LiquidFlow
    LIE
    PhaseField
    RichardsComponentTransport
    RichardsFlow
    SmallDeformation
    TES
    ThermalTwoPhaseFlowWithPP
    ThermoMechanics
    TwoPhaseFlowWithPP
    TwoPhaseFlowWithPrho)

APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES CalculateSurfaceFlux)
APPEND_SOURCE_FILES(SOURCES Output)
APPEND_SOURCE_FILES(SOURCES Parameter)
APPEND_SOURCE_FILES(SOURCES SourceTerms)
APPEND_SOURCE_FILES(SOURCES Utils)

add_library(ProcessLib ${SOURCES})

target_link_libraries(ProcessLib
    PUBLIC BaseLib GeoLib MaterialLib MathLib MeshLib NumLib logog
    PRIVATE MeshGeoToolsLib
)

if(TARGET Eigen)
    add_dependencies(ProcessLib Eigen)
endif()
if(TARGET Boost)
    add_dependencies(ProcessLib Boost)
endif()
if(OGS_INSITU)
    target_link_libraries(ProcessLib InSituLib)
endif()

if(OGS_USE_PCH)
    cotire(ProcessLib)
endif()
