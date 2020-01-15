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

target_link_libraries(
    ProcessLib
    PUBLIC BaseLib MaterialLib MathLib MeshGeoToolsLib MeshLib NumLib
    PRIVATE ParameterLib GitInfoLib spdlog::spdlog)

if (OGS_BUILD_PROCESS_ComponentTransport
    OR OGS_BUILD_PROCESS_RichardsComponentTransport)
    target_link_libraries(ProcessLib PUBLIC ChemistryLib)
endif()

if(OGS_USE_PYTHON)
    add_subdirectory(BoundaryCondition/Python)
    target_link_libraries(ProcessLib PUBLIC ProcessLibBoundaryConditionPython)

    add_subdirectory(SourceTerms/Python)
    target_link_libraries(ProcessLib PUBLIC ProcessLibSourceTermPython)
endif()

if(OGS_INSITU)
    target_link_libraries(ProcessLib PRIVATE InSituLib)
endif()
