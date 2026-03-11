add_subdirectory(SolidModels)
add_subdirectory(FractureModels)
add_subdirectory(Utils)

# Source files
append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    .
    MPL
    MPL/Properties
    MPL/Properties/CapillaryPressureSaturation
    MPL/Properties/GibbsFreeEnergy
    MPL/Properties/Density
    MPL/Properties/Enthalpy
    MPL/Properties/RelativePermeability
    MPL/Properties/SwellingStress
    MPL/Properties/ThermalConductivity
    MPL/Properties/VapourDiffusion
    MPL/Properties/Viscosity
    MPL/Components
    MPL/Utils
)

ogs_add_library(
    MaterialLib GENERATE_EXPORT_HEADER ${_sources} PUBLIC_HEADERS ${_headers}
)

target_link_libraries(
    MaterialLib
    PUBLIC BaseLib
           Eigen3::Eigen
           MaterialLib_SolidModels
           MaterialLib_FractureModels
           MaterialLib_Utils
           MathLib
           MeshLib
    PRIVATE ParameterLib exprtk Boost::math
)

target_precompile_headers(
    MaterialLib
    PRIVATE
    [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
    [["ParameterLib/Parameter.h"]]
    <Eigen/Core>
)
