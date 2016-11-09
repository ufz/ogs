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

APPEND_SOURCE_FILES(SOURCES LiquidFlow)
APPEND_SOURCE_FILES(SOURCES HydroMechanics)

APPEND_SOURCE_FILES(SOURCES TwoPhaseFlowWithPP)

APPEND_SOURCE_FILES(SOURCES SmallDeformation)

APPEND_SOURCE_FILES(SOURCES LIE/BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES LIE/Common)
APPEND_SOURCE_FILES(SOURCES LIE/HydroMechanics)
APPEND_SOURCE_FILES(SOURCES LIE/HydroMechanics/LocalAssembler)
APPEND_SOURCE_FILES(SOURCES LIE/SmallDeformation)
APPEND_SOURCE_FILES(SOURCES LIE/SmallDeformation/LocalAssembler)

APPEND_SOURCE_FILES(SOURCES TES)

add_subdirectory(HeatConduction)
APPEND_SOURCE_FILES(SOURCES HeatConduction)

APPEND_SOURCE_FILES(SOURCES RichardsFlow)

APPEND_SOURCE_FILES(SOURCES HT)

add_subdirectory(CalculateSurfaceFlux)
APPEND_SOURCE_FILES(SOURCES CalculateSurfaceFlux)

APPEND_SOURCE_FILES(SOURCES Utils)

add_library(ProcessLib ${SOURCES})

target_link_libraries(ProcessLib
    MaterialLib
    MeshGeoToolsLib
    NumLib # for shape matrices
    ${VTK_LIBRARIES}
)

ADD_VTK_DEPENDENCY(ProcessLib)

if(TARGET Eigen)
    add_dependencies(ProcessLib Eigen)
endif()
if(TARGET Boost)
    add_dependencies(ProcessLib Boost)
endif()
