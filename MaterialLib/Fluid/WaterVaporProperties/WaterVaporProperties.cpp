/**
 *  \copyright
 *   Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WaterVaporProperties.cpp
 *
 */

#include "WaterVaporProperties.h"

#include <array>
#include <cmath>

#include "MaterialLib/PhysicalConstant.h"

namespace MaterialLib
{
using PhysicalConstant::IdealGasConstant;
using PhysicalConstant::CelsiusZeroInKelvin;
namespace Fluid
{
static const double temperature_0 = 373.15;  /// reference temperature in [K]
static const double p_0 = 101325.0;          /// reference pressure
static const double h_wg = 2258000.0;  /// latent heat of water evaporation

double WaterVaporProperties::calculateSaturatedVaporPressure(
    const double T) const
{
    return p_0 * std::exp(((1 / temperature_0) - (1 / T)) * _water_mol_mass * h_wg /
                           IdealGasConstant);
}

double WaterVaporProperties::calculateDerivativedPsatdT(const double T) const
{
    return p_0 * (_water_mol_mass * h_wg / IdealGasConstant) * (1. / T / T) *
           std::exp(((1. / temperature_0) - (1. / T)) * _water_mol_mass * h_wg /
                    IdealGasConstant);
}

double WaterVaporProperties::calculateVaporPressureNonwet(
    const double pc, const double T,
    const double mass_density_water) const
{
    const double p_sat = calculateSaturatedVaporPressure(T);
    const double c_w = _water_mol_mass / IdealGasConstant / T;
    return p_sat * std::exp(-pc * c_w / mass_density_water);
}
double WaterVaporProperties::calculateDerivativedPgwdT(
    const double pc, const double T, const double mass_density_water) const
{
    const double c_w = _water_mol_mass / IdealGasConstant / T;
    const double p_sat = calculateSaturatedVaporPressure(T);
    const double dPsatdT = calculateDerivativedPsatdT(T);
    return dPsatdT * std::exp(-pc * c_w / mass_density_water) +
           p_sat * std::exp(-pc * c_w / mass_density_water) *
               (pc * _water_mol_mass / mass_density_water / IdealGasConstant / T / T);
}
double WaterVaporProperties::calculateDerivativedPgwdPC(
    const double pc, const double T, const double mass_density_water) const
{
    const double c_w = _water_mol_mass / IdealGasConstant / T;
    const double p_sat = calculateSaturatedVaporPressure(T);
    return p_sat * std::exp(-pc * c_w / mass_density_water) * (-c_w / mass_density_water);
}
double WaterVaporProperties::calculatedDensityNonwetdT(
    const double p_air_nonwet, const double p_vapor_nonwet, const double pc,
    const double T, const double mass_density_water) const
{
    const double dPgwdT = calculateDerivativedPgwdT(pc, T, mass_density_water);
    return -((p_air_nonwet * _air_mol_mass + p_vapor_nonwet * _water_mol_mass) / IdealGasConstant /
             T / T) +
           (_water_mol_mass - _air_mol_mass) * dPgwdT / IdealGasConstant / T;
}
double WaterVaporProperties::getWaterVaporEnthalpySimple(const double temperature,
    const double heat_capacity_water_vapor,
    const double /*pressure*/,
    const double /*latent_heat_evaporation*/) const
{
    return heat_capacity_water_vapor * (temperature - CelsiusZeroInKelvin) +
        h_wg;
}

}  // end namespace
}  // end namespace
