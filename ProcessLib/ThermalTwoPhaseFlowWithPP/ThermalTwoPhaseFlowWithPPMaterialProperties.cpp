/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermalTwoPhaseFlowWithPPMaterialProperties.h"

#include <logog/include/logog.hpp>

#include "BaseLib/reorderVector.h"
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
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

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
    : _two_phase_material_model(std::move(two_phase_material_model)),
      _specific_heat_capacity_solid(std::move(specific_heat_capacity_solid)),
      _specific_heat_capacity_water(std::move(specific_heat_capacity_water)),
      _specific_heat_capacity_air(std::move(specific_heat_capacity_air)),
      _specific_heat_capacity_vapor(std::move(specific_heat_capacity_vapor)),
      _thermal_conductivity_dry_solid(
          std::move(thermal_conductivity_dry_solid)),
      _thermal_conductivity_wet_solid(
          std::move(thermal_conductivity_wet_solid)),
      _water_vapor_properties(std::move(water_vapor_properties))
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
    return _specific_heat_capacity_solid->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityWater(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _specific_heat_capacity_water->getValue(vars);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityAir(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _specific_heat_capacity_air->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getSpecificHeatCapacityVapor(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _specific_heat_capacity_vapor->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getThermalConductivityDrySolid(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _thermal_conductivity_dry_solid->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getThermalConductivityWetSolid(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _thermal_conductivity_wet_solid->getValue(vars);
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateUnsatHeatConductivity(
    double const /*t*/, ProcessLib::SpatialPosition const& /*x*/,
    double const Sw, double const lambda_pm_dry,
    double const lambda_pm_wet) const
{
    double lambda_pm =
        lambda_pm_dry + std::sqrt(Sw) * (lambda_pm_wet - lambda_pm_dry);
    if (Sw > 1)
        lambda_pm = lambda_pm_wet;
    else if (Sw < 0)
        lambda_pm = lambda_pm_dry;
    return lambda_pm;
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateSaturatedVaporPressure(
    const double T) const
{
    return _water_vapor_properties->calculateSaturatedVaporPressure(T);
}
double
ThermalTwoPhaseFlowWithPPMaterialProperties::calculateVaporPressureNonwet(
    const double pc, const double T, const double mass_density_water) const
{
    return _water_vapor_properties->calculateVaporPressureNonwet(pc, T,
                                                                 mass_density_water);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPsatdT(
    const double T) const
{
    return _water_vapor_properties->calculateDerivativedPsatdT(T);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPgwdT(
    const double pc, const double T, const double mass_density_water) const
{
    return _water_vapor_properties->calculateDerivativedPgwdT(pc, T,
                                                              mass_density_water);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculateDerivativedPgwdPC(
    const double pc, const double T, const double mass_density_water) const
{
    return _water_vapor_properties->calculateDerivativedPgwdPC(pc, T,
                                                               mass_density_water);
}
double ThermalTwoPhaseFlowWithPPMaterialProperties::calculatedDensityNonwetdT(
    const double p_air_nonwet, const double p_vapor_nonwet, const double pc,
    const double T, const double mass_density_water) const
{
    return _water_vapor_properties->calculatedDensityNonwetdT(
        p_air_nonwet, p_vapor_nonwet, pc, T, mass_density_water);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getWaterVaporEnthalpySimple(
    const double temperature, const double heat_capacity_water_vapor,
    const double pg, const double latent_heat_evaporation) const
{
    return _water_vapor_properties->getWaterVaporEnthalpySimple(
        temperature, heat_capacity_water_vapor, pg, latent_heat_evaporation);
}

double ThermalTwoPhaseFlowWithPPMaterialProperties::getAirEnthalpySimple(
    const double temperature,
    const double heat_capacity_dry_air,
    const double /*pg*/) const
{
    return heat_capacity_dry_air * (temperature - CelsiusZeroInKelvin) +
           IdealGasConstant * (temperature - CelsiusZeroInKelvin) / _air_mol_mass;
}

double
ThermalTwoPhaseFlowWithPPMaterialProperties::getLiquidWaterEnthalpySimple(
    const double temperature,
    const double heat_capacity_liquid_water,
    const double /*pl*/) const
{
    return heat_capacity_liquid_water * (temperature - CelsiusZeroInKelvin);
}

}  // end of namespace
}  // end of namespace
