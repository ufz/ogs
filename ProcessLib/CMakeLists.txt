# Collect the process libraries in interface library
add_library(Processes INTERFACE)

# Check whether LIE processes are enabled
foreach(process ${_enabled_processes})
    if(process STREQUAL "LIE_M" OR process STREQUAL "LIE_HM")
        add_subdirectory(LIE/Common)
        break()
    endif()
endforeach()

foreach(process ${_enabled_processes})
    if(process STREQUAL "LIE_HM")
        add_subdirectory(LIE/HydroMechanics)
    elseif(process STREQUAL "LIE_M")
        add_subdirectory(LIE/SmallDeformation)
    else()
        add_subdirectory(${process})
    endif()
    target_link_libraries(Processes INTERFACE ${process})
    set_target_properties(${process} PROPERTIES JOB_POOL_COMPILE heavy_tasks)
endforeach()
include(Tests.cmake)

append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    .
    Assembly
    BoundaryConditionAndSourceTerm
    BoundaryConditionAndSourceTerm/Utils
    ConstitutiveRelations
    SurfaceFlux
    Output
    Utils
    Reflection
    Graph
)
# For processes that contains both hydraulical process and mechanical process
foreach(process ${_enabled_processes})
    if(process STREQUAL "LIE_HM"
       OR process STREQUAL "HydroMechanics"
       OR process STREQUAL "ThermoHydroMechanics"
       OR process STREQUAL "RichardsMechanics"
       OR process STREQUAL "TH2M"
       OR process STREQUAL "ThermoRichardsMechanics"
       OR process STREQUAL "PhaseField"
    )
        append_source_and_header_files(
            _sources _headers _private DIRECTORIES Common/HydroMechanics
        )
        break()
    endif()
endforeach()

foreach(process ${_enabled_processes})
    if(process STREQUAL "LiquidFlow" OR process STREQUAL "ThermoHydroMechanics")
        append_source_and_header_files(
            _sources _headers _private DIRECTORIES Common/HydraulicProcess
        )
        break()
    endif()
endforeach()

foreach(process ${_enabled_processes})
    if(process STREQUAL "LargeDeformation")
        append_source_and_header_files(
            _sources _headers _private DIRECTORIES Deformation
        )
        break()
    endif()
endforeach()

ogs_add_library(
    ProcessLib GENERATE_EXPORT_HEADER ${_sources} PUBLIC_HEADERS ${_headers}
)

add_subdirectory(BoundaryConditionAndSourceTerm/Python)

target_link_libraries(
    ProcessLib
    PUBLIC BaseLib
           MaterialLib
           MathLib
           MeshGeoToolsLib
           MeshLib
           NumLib
           range-v3
           $<$<TARGET_EXISTS:petsc>:petsc>
           ProcessLibBoundaryConditionAndSourceTermPython
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
