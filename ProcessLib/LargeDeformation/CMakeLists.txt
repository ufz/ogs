append_source_and_header_files(_sources _headers _private DIRECTORIES .)
append_source_and_header_files(
    _sources _headers _private DIRECTORIES ConstitutiveRelations
)

if(NOT OGS_USE_MFRONT)
    message(
        FATAL_ERROR
            "LargeDeformations process without MFront library is not available. Use enable "
            "MFront in the CMake settings via OGS_USE_MFRONT variable."
    )
endif()

ogs_add_library(LargeDeformation ${_sources} PUBLIC_HEADERS ${_headers})
target_link_libraries(LargeDeformation PUBLIC ProcessLib PRIVATE ParameterLib)

target_precompile_headers(
    LargeDeformation
    PRIVATE
    [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]]
    [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]]
    <Eigen/Core>
)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
