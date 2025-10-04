get_source_files(SOURCES)
append_source_files(SOURCES LocalAssembler)

if(WIN32)
    # Exceeds library limit of 65535 objects
    set(_static STATIC)
endif()

ogs_add_library(LIE_HM ${_static} ${SOURCES})
target_link_libraries(LIE_HM PUBLIC ProcessLib LIECommon PRIVATE ParameterLib)

target_precompile_headers(LIE_HM PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
