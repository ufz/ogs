include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Collect the process libraries in interface library
add_library(Processes INTERFACE)

foreach(process ${ProcessesList})
    if(OGS_BUILD_PROCESS_${process})
        add_subdirectory(${process})
        target_link_libraries(Processes INTERFACE ${process})
    endif()
endforeach()


APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES SurfaceFlux)
APPEND_SOURCE_FILES(SOURCES Output)
APPEND_SOURCE_FILES(SOURCES SourceTerms)
APPEND_SOURCE_FILES(SOURCES Utils)

if(MSVC)
    # Workaround for 65535 objects limit on Win DLLs, see
    # https://docs.microsoft.com/en-us/cpp/error-messages/tool-errors/linker-tools-error-lnk1189?view=vs-2017
    set(WIN_STATIC STATIC)
endif()

add_library(ProcessLib ${WIN_STATIC} ${SOURCES})
if(BUILD_SHARED_LIBS AND NOT MSVC)
    install(TARGETS ProcessLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                               RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

include(GenerateExportHeader)
generate_export_header(ProcessLib)
target_include_directories(ProcessLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(ProcessLib
    PUBLIC
        BaseLib
        MaterialLib
        MathLib
        MeshGeoToolsLib
        MeshLib
        NumLib
        logog
    PRIVATE
        ParameterLib
)

if(OGS_USE_PYTHON)
    add_subdirectory(BoundaryCondition/Python)
    target_link_libraries(ProcessLib
        PUBLIC ProcessLibBoundaryConditionPython)

    add_subdirectory(SourceTerms/Python)
    target_link_libraries(ProcessLib
        PUBLIC ProcessLibSourceTermPython)
endif()

if(OGS_INSITU)
    target_link_libraries(ProcessLib InSituLib)
endif()

if(OGS_USE_PCH)
    cotire(ProcessLib)
endif()
