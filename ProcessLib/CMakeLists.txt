# Collect the process libraries in interface library
add_library(Processes INTERFACE)

foreach(process ${_processes_list})
    if(OGS_BUILD_PROCESS_${process})
        add_subdirectory(${process})
        target_link_libraries(Processes INTERFACE ${process})
    endif()
endforeach()

get_source_files(SOURCES)
append_source_files(SOURCES BoundaryCondition)
append_source_files(SOURCES SurfaceFlux)
append_source_files(SOURCES Output)
append_source_files(SOURCES SourceTerms)
append_source_files(SOURCES Utils)

ogs_add_library(ProcessLib ${SOURCES})

if(OGS_USE_PYTHON)
    add_subdirectory(BoundaryCondition/Python)
    add_subdirectory(SourceTerms/Python)
endif()

target_link_libraries(
    ProcessLib
    PUBLIC
        BaseLib
        MaterialLib
        MathLib
        MeshGeoToolsLib
        MeshLib
        NumLib
        $<$<TARGET_EXISTS:ProcessLibBoundaryConditionPython>:ProcessLibBoundaryConditionPython>
        $<$<TARGET_EXISTS:ProcessLibSourceTermPython>:ProcessLibSourceTermPython>
        $<$<TARGET_EXISTS:petsc>:petsc>
        nlohmann_json
    PRIVATE ParameterLib GitInfoLib spdlog::spdlog
            $<$<TARGET_EXISTS:InSituLib>:InSituLib>
)

target_compile_definitions(
    ProcessLib
    PUBLIC $<$<BOOL:${OGS_USE_XDMF}>:OGS_USE_XDMF>
           # Enabled elements
           OGS_MAX_ELEMENT_DIM=${OGS_MAX_ELEMENT_DIM}
           OGS_MAX_ELEMENT_ORDER=${OGS_MAX_ELEMENT_ORDER}
           $<$<BOOL:${OGS_ENABLE_ELEMENT_SIMPLEX}>:OGS_ENABLE_ELEMENT_SIMPLEX>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_CUBOID}>:OGS_ENABLE_ELEMENT_CUBOID>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PRISM}>:OGS_ENABLE_ELEMENT_PRISM>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PYRAMID}>:OGS_ENABLE_ELEMENT_PYRAMID>
    PRIVATE $<$<BOOL:${USE_INSITU}>:USE_INSITU>
)
