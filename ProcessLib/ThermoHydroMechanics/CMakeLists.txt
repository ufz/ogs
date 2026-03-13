append_source_and_header_files(_sources _headers _private DIRECTORIES .)

ogs_add_library(
    ThermoHydroMechanics ${_sources} GENERATE_EXPORT_HEADER PUBLIC_HEADERS
    ${_headers}
)

# Disable Windows DLL export all symbols for this target to reduce symbol count
# and avoid LNK1189 error (65535 object limit exceeded)
set_target_properties(
    ThermoHydroMechanics PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS OFF
)

# Set all symbols to hidden visibility by default
set_target_properties(
    ThermoHydroMechanics
    PROPERTIES CXX_VISIBILITY_PRESET hidden C_VISIBILITY_PRESET hidden
               VISIBILITY_INLINES_HIDDEN ON
)

target_link_libraries(
    ThermoHydroMechanics PUBLIC ProcessLib PRIVATE ParameterLib
)

target_precompile_headers(
    ThermoHydroMechanics
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
