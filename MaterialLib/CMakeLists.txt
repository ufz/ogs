# Source files
get_source_files(SOURCES)
append_source_files(SOURCES Adsorption)
append_source_files(SOURCES SolidModels)
append_source_files(SOURCES FractureModels)

append_source_files(SOURCES Fluid)
append_source_files(SOURCES Fluid/Density)
append_source_files(SOURCES Fluid/Viscosity)
append_source_files(SOURCES Fluid/GibbsFreeEnergy)
append_source_files(SOURCES Fluid/FluidProperties)
append_source_files(SOURCES Fluid/SpecificHeatCapacity)
append_source_files(SOURCES Fluid/ThermalConductivity)
append_source_files(SOURCES Fluid/WaterVaporProperties)

append_source_files(SOURCES PorousMedium/Porosity)
append_source_files(SOURCES PorousMedium/Storage)
append_source_files(SOURCES PorousMedium/Permeability)
append_source_files(SOURCES PorousMedium/UnsaturatedProperty/CapillaryPressure)
append_source_files(SOURCES PorousMedium/UnsaturatedProperty/RelativePermeability)
append_source_files(SOURCES TwoPhaseModels)

add_library(MaterialLib ${SOURCES})

include(GenerateExportHeader)
generate_export_header(MaterialLib)
target_include_directories(MaterialLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(MaterialLib
    PUBLIC BaseLib NumLib logog # ProcessLib
    PRIVATE MathLib MeshLib
)

if(OGS_USE_PCH)
    cotire(MaterialLib)
endif()
