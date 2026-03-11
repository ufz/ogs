append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    .
    BHE
    BoundaryConditions
    LocalAssemblers
)

ogs_add_library(HeatTransportBHE ${_sources} PUBLIC_HEADERS ${_headers})

target_link_libraries(
    HeatTransportBHE PUBLIC ProcessLib
    PRIVATE ParameterLib
            $<$<TARGET_EXISTS:pybind11::pybind11>:pybind11::module>
)

target_precompile_headers(HeatTransportBHE PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
