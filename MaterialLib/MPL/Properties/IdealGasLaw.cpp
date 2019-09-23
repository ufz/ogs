/**
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/IdealGasLaw.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
double molarMass(Phase* _phase, Component* _component,
                 VariableArray const& variable_array,
                 ParameterLib::SpatialPosition const& pos, double const t)
{
    if (_phase)  // IdealGasLaw of an entire phase
    {
        return _phase->property(PropertyType::molar_mass)
            .template value<double>(variable_array, pos, t);
    }
    if (_component)  // IdealGasLaw of a single component
    {
        return _component->property(PropertyType::molar_mass)
            .template value<double>(variable_array, pos, t);
    }
    OGS_FATAL(
        "Neither a phase nor a component are set for retrieving molar "
        "mass.");

    return 0.;
}

PropertyDataType IdealGasLaw::value(VariableArray const& variable_array,
                                    ParameterLib::SpatialPosition const& pos,
                                    double const t) const
{
    const double gas_constant = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pressure = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    const double temperature = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    double molar_mass = molarMass(_phase, _component, variable_array, pos, t);

    const double density = pressure * molar_mass / gas_constant / temperature;

    return density;
}

PropertyDataType IdealGasLaw::dValue(VariableArray const& variable_array,
                                     Variable const primary_variable,
                                     ParameterLib::SpatialPosition const& pos,
                                     double const t) const
{
    const double gas_constant = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pressure = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    const double temperature = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    double molar_mass = molarMass(_phase, _component, variable_array, pos, t);

    if (primary_variable == Variable::temperature)
    {
        return -pressure * molar_mass / gas_constant / temperature /
               temperature;
    }

    if (primary_variable == Variable::phase_pressure)
    {
        return molar_mass / gas_constant / temperature;
    }

    OGS_FATAL(
        "IdealGasLaw::dValue is implemented for derivatives "
        "with respect to phase pressure or temperature only.");

    return 0.;
}

PropertyDataType IdealGasLaw::d2Value(VariableArray const& variable_array,
                                      Variable const primary_variable1,
                                      Variable const primary_variable2,
                                      ParameterLib::SpatialPosition const& pos,
                                      double const t) const
{
    const double gas_constant = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pressure = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    const double temperature = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    double molar_mass = molarMass(_phase, _component, variable_array, pos, t);

    if ((primary_variable1 == Variable::phase_pressure) &&
        (primary_variable2 == Variable::phase_pressure))
    {
        // d2rho_dp2
        return 0.;
    }
    if ((primary_variable1 == Variable::temperature) &&
        (primary_variable2 == Variable::temperature))
    {
        // d2rho_dT2
        return 2. * molar_mass * pressure / gas_constant / temperature /
               temperature / temperature;
    }
    if (((primary_variable1 == Variable::phase_pressure) &&
         (primary_variable2 == Variable::temperature)) ||
        ((primary_variable1 == Variable::temperature) &&
         (primary_variable2 == Variable::phase_pressure)))
    {
        // d2rho_dpdT or d2rho_dTdp
        return -molar_mass / gas_constant / temperature / temperature;
    }

    OGS_FATAL(
        "IdealGasLaw::d2Value is implemented for derivatives "
        "with respect to phase pressure and temperature only.");

    return 0.;
}

}  // namespace MaterialPropertyLib
