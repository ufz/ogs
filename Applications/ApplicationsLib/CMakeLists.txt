# Source files
get_source_files(SOURCES_APPLICATIONSLIB)
set(LIB_SOURCES ${SOURCES_APPLICATIONSLIB})

# Library
ogs_add_library(ApplicationsLib ${LIB_SOURCES})

target_link_libraries(
    ApplicationsLib
    PUBLIC BaseLib GeoLib NumLib Processes
    PRIVATE CMakeInfoLib
            MathLib
            MeshLib
            MeshGeoToolsLib
            ParameterLib
            ProcessLib
            spdlog::spdlog
            $<$<BOOL:${OGS_USE_PYTHON}>:pybind11::pybind11>
            $<$<BOOL:${OGS_USE_PETSC}>:petsc>
)

if(OGS_USE_LIS)
    target_include_directories(ApplicationsLib PUBLIC ${LIS_INCLUDE_DIR})
endif()

target_compile_definitions(ApplicationsLib PUBLIC $<$<BOOL:${USE_LIS}>:USE_LIS>)

# Set cpp definitions if the cmake option is enabled for the given process.
foreach(process ${_processes_list})
    if(OGS_BUILD_PROCESS_${process})
        string(TOUPPER "OGS_BUILD_PROCESS_${process}" EnableProcess)
        set_property(
            TARGET ApplicationsLib APPEND PROPERTY COMPILE_DEFINITIONS
                                                   ${EnableProcess}
        )
    endif()
endforeach()
