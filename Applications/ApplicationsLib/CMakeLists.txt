# Source files
get_source_files(SOURCES_APPLICATIONSLIB)
set(LIB_SOURCES ${SOURCES_APPLICATIONSLIB})

# Library
ogs_add_library(ApplicationsLib ${LIB_SOURCES})

target_link_libraries(
    ApplicationsLib
    PUBLIC BaseLib GeoLib NumLib
           $<$<TARGET_EXISTS:VTK::ParallelMPI>:VTK::ParallelMPI>
    PRIVATE CMakeInfoLib
            MathLib
            MeshLib
            MeshGeoToolsLib
            MeshToolsLib
            pybind11::pybind11
            $<$<BOOL:${OGS_USE_PETSC}>:PkgConfig::PETSC>
)

if(OGS_BUILD_CLI OR OGS_BUILD_UTILS OR OGS_BUILD_TESTING)
    target_link_libraries(
        ApplicationsLib PUBLIC Processes PRIVATE ParameterLib ProcessLib
    )
elseif(OGS_BUILD_GUI)
    target_link_libraries(
        ApplicationsLib PRIVATE MeshToolsLib nlohmann_json::nlohmann_json
    )
endif()

if(OGS_USE_LIS)
    target_include_directories(ApplicationsLib PUBLIC ${LIS_INCLUDE_DIR})
endif()

target_compile_definitions(
    ApplicationsLib PUBLIC $<$<BOOL:${USE_LIS}>:USE_LIS>
    PRIVATE $<$<BOOL:${OGS_USE_INSITU}>:OGS_USE_INSITU>
)

# ~~~
# Prevent warning:
# warning: ‘ProjectData’ declared with greater visibility than the type of its field ‘ProjectData::_py_scoped_interpreter’
# see https://pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
# ~~~
if(NOT OGS_BUILD_WHEEL AND OGS_BUILD_CLI)
    target_link_libraries(ApplicationsLib PRIVATE ogs_embedded_python)
    set_target_properties(
        ApplicationsLib PROPERTIES CXX_VISIBILITY_PRESET hidden
    )
endif()

# Set cpp definitions if the cmake option is enabled for the given process.
foreach(process ${_enabled_processes})
    string(TOUPPER "OGS_BUILD_PROCESS_${process}" EnableProcess)
    set_property(
        TARGET ApplicationsLib APPEND PROPERTY COMPILE_DEFINITIONS
                                               ${EnableProcess}
    )
endforeach()

target_precompile_headers(
    ApplicationsLib PRIVATE [["BaseLib/Error.h"]] [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
)
