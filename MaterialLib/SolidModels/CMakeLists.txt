get_source_files(SOURCES)

ogs_add_library(MaterialLib_SolidModels ${SOURCES})
target_link_libraries(
    MaterialLib_SolidModels PUBLIC MathLib Boost::math
    PRIVATE ParameterLib MaterialLib_Utils NumLib
)

if(OGS_USE_MFRONT)
    add_subdirectory(MFront)
    target_compile_definitions(MaterialLib_SolidModels PUBLIC OGS_USE_MFRONT)
    target_link_libraries(
        MaterialLib_SolidModels PUBLIC MaterialLib_SolidModels_MFront
    )
endif()

target_precompile_headers(
    MaterialLib_SolidModels PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]] <Eigen/Core>
)
