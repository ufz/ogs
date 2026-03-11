append_source_and_header_files(_sources _headers _private DIRECTORIES .)

if(WIN32)
    # Exceeds library limit of 65535 objects
    set(_static STATIC)
endif()

ogs_add_library(
    ComponentTransport ${_static} ${_sources} PUBLIC_HEADERS ${_headers}
)

target_link_libraries(
    ComponentTransport PUBLIC ProcessLib ChemistryLib PRIVATE ParameterLib
)

target_precompile_headers(ComponentTransport PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
