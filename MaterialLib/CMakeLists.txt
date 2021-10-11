add_subdirectory(SolidModels)
add_subdirectory(FractureModels)

# Source files
get_source_files(SOURCES)
append_source_files(SOURCES Adsorption)

append_source_files(SOURCES Fluid)
append_source_files(SOURCES Fluid/Density)
append_source_files(SOURCES Fluid/Viscosity)
append_source_files(SOURCES Fluid/GibbsFreeEnergy)
append_source_files(SOURCES Fluid/FluidProperties)
append_source_files(SOURCES Fluid/SpecificHeatCapacity)
append_source_files(SOURCES Fluid/ThermalConductivity)
append_source_files(SOURCES Fluid/WaterVaporProperties)

append_source_files(SOURCES MPL)
append_source_files(SOURCES MPL/Properties)
append_source_files(SOURCES MPL/Properties/CapillaryPressureSaturation)
append_source_files(SOURCES MPL/Properties/Density)
append_source_files(SOURCES MPL/Properties/Enthalpy)
append_source_files(SOURCES MPL/Properties/RelativePermeability)
append_source_files(SOURCES MPL/Properties/SwellingStress)
append_source_files(SOURCES MPL/Properties/ThermalConductivity)
append_source_files(SOURCES MPL/Properties/VapourDiffusion)
append_source_files(SOURCES MPL/Components)
append_source_files(SOURCES MPL/Utils)

append_source_files(SOURCES PorousMedium)
append_source_files(SOURCES PorousMedium/Porosity)
append_source_files(SOURCES PorousMedium/Storage)
append_source_files(SOURCES PorousMedium/Permeability)
append_source_files(SOURCES PorousMedium/UnsaturatedProperty/CapillaryPressure)
append_source_files(
    SOURCES PorousMedium/UnsaturatedProperty/RelativePermeability
)
append_source_files(SOURCES TwoPhaseModels)

ogs_add_library(MaterialLib GENERATE_EXPORT_HEADER ${SOURCES})

target_link_libraries(
    MaterialLib PUBLIC MaterialLib_SolidModels MaterialLib_FractureModels
    PRIVATE MathLib MeshLib ParameterLib exprtk
)
