include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Source files
APPEND_SOURCE_FILES(SOURCES)

add_subdirectory(BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES BoundaryCondition)

APPEND_SOURCE_FILES(SOURCES Deformation)

add_subdirectory(Parameter)
APPEND_SOURCE_FILES(SOURCES Parameter)

add_subdirectory(GroundwaterFlow)
APPEND_SOURCE_FILES(SOURCES GroundwaterFlow)

add_subdirectory(LiquidFlow)
APPEND_SOURCE_FILES(SOURCES LiquidFlow)
add_subdirectory(HydroMechanics)
APPEND_SOURCE_FILES(SOURCES HydroMechanics)

add_subdirectory(TwoPhaseFlowWithPP)
APPEND_SOURCE_FILES(SOURCES TwoPhaseFlowWithPP)
add_subdirectory(TwoPhaseFlowWithPrho)
APPEND_SOURCE_FILES(SOURCES TwoPhaseFlowWithPrho)

add_subdirectory(SmallDeformation)
APPEND_SOURCE_FILES(SOURCES SmallDeformation)

add_subdirectory(LIE/SmallDeformation)
add_subdirectory(LIE/HydroMechanics)
APPEND_SOURCE_FILES(SOURCES LIE/BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES LIE/Common)
APPEND_SOURCE_FILES(SOURCES LIE/HydroMechanics)
APPEND_SOURCE_FILES(SOURCES LIE/HydroMechanics/LocalAssembler)
APPEND_SOURCE_FILES(SOURCES LIE/SmallDeformation)
APPEND_SOURCE_FILES(SOURCES LIE/SmallDeformation/LocalAssembler)

add_subdirectory(TES)
APPEND_SOURCE_FILES(SOURCES TES)

add_subdirectory(HeatConduction)
APPEND_SOURCE_FILES(SOURCES HeatConduction)

add_subdirectory(RichardsFlow)
APPEND_SOURCE_FILES(SOURCES RichardsFlow)

add_subdirectory(HC)
APPEND_SOURCE_FILES(SOURCES HC)

add_subdirectory(HT)
APPEND_SOURCE_FILES(SOURCES HT)

add_subdirectory(CalculateSurfaceFlux)
APPEND_SOURCE_FILES(SOURCES CalculateSurfaceFlux)

APPEND_SOURCE_FILES(SOURCES Utils)

add_library(ProcessLib ${SOURCES})

target_link_libraries(ProcessLib
    MaterialLib
    MeshGeoToolsLib
    MeshLib
    NumLib # for shape matrices
)

ADD_VTK_DEPENDENCY(ProcessLib)

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
