# Source files
get_source_files(SOURCES)
append_source_files(SOURCES Adsorption)
append_source_files(SOURCES SolidModels)
append_source_files(SOURCES FractureModels)
append_source_files(SOURCES FractureModels/Permeability)

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

if(OGS_USE_MFRONT)
    add_subdirectory(SolidModels/MFront)
endif()

ogs_add_library(MaterialLib ${SOURCES})

if(OGS_USE_MFRONT)
    target_compile_definitions(MaterialLib PUBLIC OGS_USE_MFRONT)
endif()
target_link_libraries(
    MaterialLib
    PUBLIC
        $<$<TARGET_EXISTS:MaterialLib_SolidModels_MFront>:MaterialLib_SolidModels_MFront>
    PRIVATE MathLib MeshLib ParameterLib spdlog::spdlog
)
