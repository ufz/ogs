get_source_files(SOURCES)
append_source_files(SOURCES ConstitutiveRelations)

if(NOT OGS_USE_MFRONT)
    message(
        FATAL_ERROR
            "LargeDeformations process without MFront library is not available. Use enable "
            "MFront in the CMake settings via OGS_USE_MFRONT variable."
    )
endif()

ogs_add_library(LargeDeformation ${SOURCES})
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
