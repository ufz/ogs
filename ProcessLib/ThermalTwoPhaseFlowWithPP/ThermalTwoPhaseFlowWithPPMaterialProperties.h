/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/Fluid/WaterVaporProperties/WaterVaporProperties.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MaterialLib/TwoPhaseModels/TwoPhaseFlowWithPPMaterialProperties.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}
/** TODO
* in this implementation, the thermal properties of gas component directly use
* the properties of air,
* e.g. the specific heat capacity of gas component use the specific heat
* capacity of air, and since a constant specific heat capacity of air is
* assumed, the enthalpy of air is calculated based on a simplified model.
* Next, a strategy which can automatically(or from input file)switch between
* different gas components is required, and should be implemented as a material
* class, also different enthalpy models for different components are need to be
* implemented.
*/
namespace ProcessLib
{
class SpatialPosition;
namespace ThermalTwoPhaseFlowWithPP
{
class ThermalTwoPhaseFlowWithPPMaterialProperties final
{
public:
    using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

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
            water_vapor_properties);

    double getSpecificHeatCapacitySolid(const double p, const double T) const;
    double getSpecificHeatCapacityWater(const double p, const double T) const;
    double getSpecificHeatCapacityAir(const double p, const double T) const;
    double getSpecificHeatCapacityVapor(const double p, const double T) const;
    double getThermalConductivityDrySolid(const double p, const double T) const;
    double getThermalConductivityWetSolid(const double p, const double T) const;
    /// Calculates the unsaturated heat conductivity
    double calculateUnsatHeatConductivity(double const t,
                                          ProcessLib::SpatialPosition const& x,
                                          double const Sw,
                                          double const lambda_pm_dry,
                                          double const lambda_pm_wet) const;
    /// water vapor saturation pressure
    double calculateSaturatedVaporPressure(const double T) const;
    /// partial water vapor pressure in nonwetting phase
    /// Kelvin equation
    double calculateVaporPressureNonwet(const double pc, const double T,
                                        const double mass_density_water) const;
    /// Derivative of SaturatedVaporPressure in terms of T
    double calculateDerivativedPsatdT(const double T) const;
    /// Derivative of partial vapor pressure in terms of T
    double calculateDerivativedPgwdT(const double pc, const double T,
                                     const double mass_density_water) const;
    /// Derivative of partial vapor pressure in terms of PC
    double calculateDerivativedPgwdPC(const double pc, const double T,
                                      const double mass_density_water) const;
    ///
    double calculatedDensityNonwetdT(const double p_air_nonwet,
                                 const double p_vapor_nonwetconst, double pc,
                                 const double T,
                                 const double mass_density_water) const;
    /// Specific enthalpy of water vapor
    double getWaterVaporEnthalpySimple(
        const double temperature,
        const double heat_capacity_water_vapor,
        const double pg,
        const double latent_heat_evaporation) const;
    /// Specific enthalpy of air
    double getAirEnthalpySimple(const double temperature,
                                const double heat_capacity_water_air,
                                const double /*pg*/) const;
    /// Specific enthalpy of liquid water
    double getLiquidWaterEnthalpySimple(const double temperature,
                                        const double heat_capacity_liquid_water,
                                        const double /*pl*/) const;
    const MaterialLib::TwoPhaseFlowWithPP::TwoPhaseFlowWithPPMaterialProperties&
    getTwoPhaseMaterialModel() const
    {
        return *_two_phase_material_model;
    }

private:
    double const& _air_mol_mass = MaterialLib::PhysicalConstant::MolarMass::Air;
    std::unique_ptr<MaterialLib::TwoPhaseFlowWithPP::
                        TwoPhaseFlowWithPPMaterialProperties> const
        _two_phase_material_model;

    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const
        _specific_heat_capacity_solid;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const
        _specific_heat_capacity_water;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const
        _specific_heat_capacity_air;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const
        _specific_heat_capacity_vapor;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const
        _thermal_conductivity_dry_solid;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const
        _thermal_conductivity_wet_solid;
    std::unique_ptr<MaterialLib::Fluid::WaterVaporProperties> const
        _water_vapor_properties;
};

}  // end of namespace
}  // end of namespace
