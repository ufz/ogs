append_source_and_header_files(_sources _headers _private DIRECTORIES .)
append_source_and_header_files(
    _sources _headers _private DIRECTORIES ConstitutiveRelations
)

ogs_add_library(SmallDeformation ${_sources} PUBLIC_HEADERS ${_headers})
target_link_libraries(SmallDeformation PUBLIC ProcessLib PRIVATE ParameterLib)

target_precompile_headers(SmallDeformation PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
