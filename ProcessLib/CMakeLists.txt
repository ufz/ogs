# Collect the process libraries in interface library
add_library(Processes INTERFACE)

foreach(process ${_enabled_processes})
    add_subdirectory(${process})
    target_link_libraries(Processes INTERFACE ${process})
    set_target_properties(${process} PROPERTIES JOB_POOL_COMPILE heavy_tasks)
endforeach()
include(Tests.cmake)

get_source_files(SOURCES)
append_source_files(SOURCES Assembly)
append_source_files(SOURCES BoundaryConditionAndSourceTerm)
append_source_files(SOURCES BoundaryConditionAndSourceTerm/Utils)
append_source_files(SOURCES ConstitutiveRelations)
append_source_files(SOURCES SurfaceFlux)
append_source_files(SOURCES Output)
append_source_files(SOURCES Utils)
append_source_files(SOURCES Reflection)
append_source_files(SOURCES Graph)
### For processes that contains both hydraulical process and mechanical process
# Define a function to check if enabled_processes contains LIE, HydroMechanics,
# ThermoHydroMechanics, RichardsMechanics, or ThermoRichardsMechanics
function(check_enabled_processes contains_hm_process)
    set(${contains_hm_process} FALSE PARENT_SCOPE)
    foreach(process ${_enabled_processes})
        if(process STREQUAL "LIE"
             OR process STREQUAL "HydroMechanics"
             OR process STREQUAL "ThermoHydroMechanics"
             OR process STREQUAL "RichardsMechanics"
             OR process STREQUAL "TH2M"
             OR process STREQUAL "ThermoRichardsMechanics")
            set(${contains_hm_process} TRUE PARENT_SCOPE)
            return()
        endif()
    endforeach()
endfunction()
check_enabled_processes(contains_hm_process)
if(contains_hm_process)
    append_source_files(SOURCES Common/HydroMechanics)
endif()

ogs_add_library(ProcessLib GENERATE_EXPORT_HEADER ${SOURCES})

add_subdirectory(BoundaryConditionAndSourceTerm/Python)

target_link_libraries(
    ProcessLib
    PUBLIC
        BaseLib
        MaterialLib
        MathLib
        MeshGeoToolsLib
        MeshLib
        NumLib
        Boost::boost
        range-v3
        $<$<TARGET_EXISTS:ProcessLibBoundaryConditionAndSourceTermPython>:ProcessLibBoundaryConditionAndSourceTermPython>
        $<$<TARGET_EXISTS:petsc>:petsc>
    PRIVATE ParameterLib GitInfoLib $<$<TARGET_EXISTS:InSituLib>:InSituLib>
)

target_compile_definitions(
    ProcessLib
    PUBLIC # Enabled elements
           OGS_MAX_ELEMENT_DIM=${OGS_MAX_ELEMENT_DIM}
           OGS_MAX_ELEMENT_ORDER=${OGS_MAX_ELEMENT_ORDER}
           $<$<BOOL:${OGS_ENABLE_ELEMENT_SIMPLEX}>:OGS_ENABLE_ELEMENT_SIMPLEX>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_CUBOID}>:OGS_ENABLE_ELEMENT_CUBOID>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PRISM}>:OGS_ENABLE_ELEMENT_PRISM>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PYRAMID}>:OGS_ENABLE_ELEMENT_PYRAMID>
    PRIVATE $<$<BOOL:${OGS_USE_INSITU}>:OGS_USE_INSITU>
)

target_precompile_headers(
    ProcessLib
    PRIVATE
    [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
    [["ParameterLib/Parameter.h"]]
    <Eigen/Core>
)
