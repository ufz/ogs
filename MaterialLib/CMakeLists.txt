add_subdirectory(SolidModels)
add_subdirectory(FractureModels)
add_subdirectory(Utils)

# Source files
get_source_files(SOURCES)

append_source_files(SOURCES MPL)
append_source_files(SOURCES MPL/Properties)
append_source_files(SOURCES MPL/Properties/CapillaryPressureSaturation)
append_source_files(SOURCES MPL/Properties/GibbsFreeEnergy)
append_source_files(SOURCES MPL/Properties/Density)
append_source_files(SOURCES MPL/Properties/Enthalpy)
append_source_files(SOURCES MPL/Properties/RelativePermeability)
append_source_files(SOURCES MPL/Properties/SwellingStress)
append_source_files(SOURCES MPL/Properties/ThermalConductivity)
append_source_files(SOURCES MPL/Properties/VapourDiffusion)
append_source_files(SOURCES MPL/Properties/Viscosity)
append_source_files(SOURCES MPL/Components)
append_source_files(SOURCES MPL/Utils)

append_source_files(SOURCES PorousMedium/Permeability)
append_source_files(SOURCES PorousMedium/UnsaturatedProperty/CapillaryPressure)
append_source_files(
    SOURCES PorousMedium/UnsaturatedProperty/RelativePermeability
)

ogs_add_library(MaterialLib GENERATE_EXPORT_HEADER ${SOURCES})

target_link_libraries(
    MaterialLib PUBLIC BaseLib Eigen3::Eigen MaterialLib_SolidModels
                       MaterialLib_FractureModels MaterialLib_Utils
    PRIVATE MathLib MeshLib ParameterLib exprtk Boost::math
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
