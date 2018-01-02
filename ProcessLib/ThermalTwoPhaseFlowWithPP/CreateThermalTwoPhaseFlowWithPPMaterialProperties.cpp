/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermalTwoPhaseFlowWithPPMaterialProperties.h"

#include <logog/include/logog.hpp>
#include <tuple>

#include "BaseLib/reorderVector.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/SpecificHeatCapacity/CreateSpecificFluidHeatCapacityModel.h"
#include "MaterialLib/Fluid/ThermalConductivity/CreateFluidThermalConductivityModel.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MaterialLib/TwoPhaseModels/CreateTwoPhaseFlowMaterialProperties.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

#include "ThermalTwoPhaseFlowWithPPMaterialProperties.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties>
createThermalTwoPhaseFlowWithPPMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const& material_ids,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    DBUG(
        "Reading material properties of nonisothermal two-phase flow process.");
    auto two_phase_model_tuple =
        MaterialLib::TwoPhaseFlowWithPP::createTwoPhaseFlowMaterialProperties(
            config, material_ids, parameters);
    auto two_phase_material_model =
        std::move(std::get<0>(two_phase_model_tuple));
    auto const& fluid_config = std::get<1>(two_phase_model_tuple);

    // Get fluid properties
    auto const& spec_heat_capacity_solid_conf =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property__specific_heat_capacity_solid}
        fluid_config.getConfigSubtree("specific_heat_capacity_solid");
    auto specific_heat_capacity_solid =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(
            spec_heat_capacity_solid_conf);
    auto const& spec_heat_capacity_water_conf =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property__specific_heat_capacity_water}
        fluid_config.getConfigSubtree("specific_heat_capacity_water");
    auto specific_heat_capacity_water =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(
            spec_heat_capacity_water_conf);
    auto const& spec_heat_capacity_air_conf =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property__specific_heat_capacity_air}
        fluid_config.getConfigSubtree("specific_heat_capacity_air");
    auto specific_heat_capacity_air =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(
            spec_heat_capacity_air_conf);
    auto const& spec_heat_capacity_vapor_conf =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property__specific_heat_capacity_water_vapor}
        fluid_config.getConfigSubtree("specific_heat_capacity_water_vapor");
    auto specific_heat_capacity_vapor =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(
            spec_heat_capacity_vapor_conf);

    auto const& thermal_conductivity_dry_solid_conf =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property__thermal_conductivity_dry_solid}
        fluid_config.getConfigSubtree("thermal_conductivity_dry_solid");
    auto thermal_conductivity_dry_solid =
        MaterialLib::Fluid::createFluidThermalConductivityModel(
            thermal_conductivity_dry_solid_conf);
    auto const& thermal_conductivity_wet_solid_conf =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property__thermal_conductivity_wet_solid}
        fluid_config.getConfigSubtree("thermal_conductivity_wet_solid");
    auto thermal_conductivity_wet_solid =
        MaterialLib::Fluid::createFluidThermalConductivityModel(
            thermal_conductivity_wet_solid_conf);

    std::unique_ptr<MaterialLib::Fluid::WaterVaporProperties> vapor_property =
        std::make_unique<MaterialLib::Fluid::WaterVaporProperties>();

    return std::make_unique<ThermalTwoPhaseFlowWithPPMaterialProperties>(
        std::move(two_phase_material_model),
        std::move(specific_heat_capacity_solid),
        std::move(specific_heat_capacity_water),
        std::move(specific_heat_capacity_air),
        std::move(specific_heat_capacity_vapor),
        std::move(thermal_conductivity_dry_solid),
        std::move(thermal_conductivity_wet_solid),
        std::move(vapor_property));
}

}  // end namespace
}  // end namespace
