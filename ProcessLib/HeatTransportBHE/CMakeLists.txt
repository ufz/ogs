get_source_files(SOURCES)
append_source_files(SOURCES BHE)
append_source_files(SOURCES BoundaryConditions)
append_source_files(SOURCES LocalAssemblers)

ogs_add_library(HeatTransportBHE ${SOURCES})

target_link_libraries(
    HeatTransportBHE PUBLIC ProcessLib
    PRIVATE ParameterLib
            $<$<TARGET_EXISTS:pybind11::pybind11>:pybind11::pybind11>
)

target_precompile_headers(HeatTransportBHE PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
