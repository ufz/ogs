include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Collect the process libraries in interface library
add_library(Processes INTERFACE)

foreach(process ${ProcessesList})
    if(OGS_BUILD_PROCESS_${process})
        add_subdirectory(${process})
        target_link_libraries(Processes INTERFACE ${process})
    endif()
endforeach()

append_source_files(SOURCES)
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
if(OGS_USE_PETSC)
    target_compile_definitions(ProcessLib PUBLIC USE_PETSC)
endif()
if(OGS_USE_XDMF)
    target_compile_definitions(ProcessLib PUBLIC OGS_USE_XDMF)
endif()
target_link_libraries(ProcessLib
    PUBLIC
        BaseLib
        MaterialLib
        MathLib
        MeshGeoToolsLib
        MeshLib
        NumLib
        $<$<OR:$<BOOL:${OGS_BUILD_PROCESS_ComponentTransport}>,$<BOOL:${OGS_BUILD_PROCESS_RichardsComponentTransport}>>:ChemistryLib>
        $<$<TARGET_EXISTS:ProcessLibBoundaryConditionPython>:ProcessLibBoundaryConditionPython>
        $<$<TARGET_EXISTS:ProcessLibSourceTermPython>:ProcessLibSourceTermPython>
    PRIVATE
        ParameterLib
        GitInfoLib
        spdlog::spdlog
        $<$<TARGET_EXISTS:InSituLib>:InSituLib>
)
