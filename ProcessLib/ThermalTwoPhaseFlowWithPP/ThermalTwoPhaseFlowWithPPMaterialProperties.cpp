/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermalTwoPhaseFlowWithPPMaterialProperties.h"

#include "BaseLib/Logging.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MaterialLib/TwoPhaseModels/TwoPhaseFlowWithPPMaterialProperties.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib
{
using MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
using MaterialLib::PhysicalConstant::IdealGasConstant;

namespace ThermalTwoPhaseFlowWithPP
{
ThermalTwoPhaseFlowWithPPMaterialProperties::
    ThermalTwoPhaseFlowWithPPMaterialProperties(
        std::unique_ptr<MaterialLib::TwoPhaseFlowWithPP::
                            TwoPhaseFlowWithPPMaterialProperties>&&
            two_phase_material_model,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            specific_heat_capacity_solid,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            specific_heat_capacity_water,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            specific_heat_capacity_air,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            specific_heat_capacity_vapor,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            thermal_conductivity_dry_solid,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            thermal_conductivity_wet_solid,
        std::unique_ptr<MaterialLib::Fluid::WaterVaporProperties>&&
            water_vapor_properties)
    : two_phase_material_model_(std::move(two_phase_material_model)),
      specific_heat_capacity_solid_(std::move(specific_heat_capacity_solid)),
      specific_heat_capacity_water_(std::move(specific_heat_capacity_water)),
      specific_heat_capacity_air_(std::move(specific_heat_capacity_air)),
      specific_heat_capacity_vapor_(std::move(specific_heat_capacity_vapor)),
      thermal_conductivity_dry_solid_(
          std::move(thermal_conductivity_dry_solid)),
      thermal_conductivity_wet_solid_(
          std::move(thermal_conductivity_wet_solid)),
      water_vapor_properties_(std::move(water_vapor_properties))
{
    DBUG("Create material properties for non-isothermal two-phase flow model.");
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacitySolid(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return specific_heat_capacity_solid_->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityWater(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return specific_heat_capacity_water_->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityAir(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return specific_heat_capacity_air_->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityVapor(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return specific_heat_capacity_vapor_->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getThermalConductivityDrySolid(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return thermal_conductivity_dry_solid_->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getThermalConductivityWetSolid(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return thermal_conductivity_wet_solid_->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateUnsatHeatConductivity(
    double const /*t*/, ParameterLib::SpatialPosition const& /*x*/,
    double const Sw, double const lambda_pm_dry,
    double const lambda_pm_wet) const
{
    double lambda_pm =
        lambda_pm_dry + std::sqrt(Sw) * (lambda_pm_wet - lambda_pm_dry);
    if (Sw > 1)
    {
        lambda_pm = lambda_pm_wet;
    }
    else if (Sw < 0)
    {
        lambda_pm = lambda_pm_dry;
    }
    return lambda_pm;
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateSaturatedVaporPressure(
    const double T) const
{
    return water_vapor_properties_->calculateSaturatedVaporPressure(T);
}
double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateVaporPressureNonwet(
    const double pc, const double T, const double mass_density_water) const
{
    return water_vapor_properties_->calculateVaporPressureNonwet(pc, T,
                                                                 mass_density_water);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPsatdT(
    const double T) const
{
    return water_vapor_properties_->calculateDerivativedPsatdT(T);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPgwdT(
    const double pc, const double T, const double mass_density_water) const
{
    return water_vapor_properties_->calculateDerivativedPgwdT(pc, T,
                                                              mass_density_water);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPgwdPC(
    const double pc, const double T, const double mass_density_water) const
{
    return water_vapor_properties_->calculateDerivativedPgwdPC(pc, T,
                                                               mass_density_water);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculatedDensityNonwetdT(
    const double p_air_nonwet, const double p_vapor_nonwet, const double pc,
    const double T, const double mass_density_water) const
{
    return water_vapor_properties_->calculatedDensityNonwetdT(
        p_air_nonwet, p_vapor_nonwet, pc, T, mass_density_water);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getWaterVaporEnthalpySimple(
    const double temperature, const double heat_capacity_water_vapor,
    const double pg, const double latent_heat_evaporation) const
{
    return water_vapor_properties_->getWaterVaporEnthalpySimple(
        temperature, heat_capacity_water_vapor, pg, latent_heat_evaporation);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getAirEnthalpySimple(
    const double temperature,
    const double heat_capacity_dry_air,
    const double /*pg*/) const
{
    return heat_capacity_dry_air * (temperature - CelsiusZeroInKelvin) +
           IdealGasConstant * (temperature - CelsiusZeroInKelvin) / air_mol_mass_;
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getLiquidWaterEnthalpySimple(
    const double temperature,
    const double heat_capacity_liquid_water,
    const double /*pl*/) const
{
    return heat_capacity_liquid_water * (temperature - CelsiusZeroInKelvin);
}

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
