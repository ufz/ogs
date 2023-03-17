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
append_source_files(SOURCES SurfaceFlux)
append_source_files(SOURCES Output)
append_source_files(SOURCES Utils)
append_source_files(SOURCES Reflection)

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
