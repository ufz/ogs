/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
IdealGasLaw::IdealGasLaw(std::string name)
{
    name_ = std::move(name);
}

PropertyDataType IdealGasLaw::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double gas_constant = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pressure = variable_array.gas_phase_pressure;
    const double temperature = variable_array.temperature;
    const double molar_mass = variable_array.molar_mass;

    const double density = pressure * molar_mass / gas_constant / temperature;

    return density;
}

PropertyDataType IdealGasLaw::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double gas_constant = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pressure = variable_array.gas_phase_pressure;
    const double temperature = variable_array.temperature;
    const double molar_mass = variable_array.molar_mass;
    // todo: add molar mass derivatives

    if (variable == Variable::temperature)
    {
        // extend to take temperature-dependent molar mass into account
        return -pressure * molar_mass / gas_constant / temperature /
               temperature;
    }

    if (variable == Variable::gas_phase_pressure)
    {
        // extend to take pressure-dependent molar mass into account
        return molar_mass / gas_constant / temperature;
    }

    OGS_FATAL(
        "IdealGasLaw::dValue is implemented for derivatives with respect to "
        "phase pressure or temperature only.");

    return 0.;
}

PropertyDataType IdealGasLaw::d2Value(
    VariableArray const& variable_array, Variable const variable1,
    Variable const variable2, ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/, double const /*dt*/) const
{
    const double gas_constant = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pressure = variable_array.gas_phase_pressure;
    const double temperature = variable_array.temperature;
    const double molar_mass = variable_array.molar_mass;
    // todo: add molar mass derivatives

    if ((variable1 == Variable::gas_phase_pressure) &&
        (variable2 == Variable::gas_phase_pressure))
    {
        // d2rho_dp2
        // extend to take pressure-dependent molar mass into account
        return 0.;
    }
    if ((variable1 == Variable::temperature) &&
        (variable2 == Variable::temperature))
    {
        // d2rho_dT2
        // extend to take temperature-dependent molar mass into account
        return 2. * molar_mass * pressure / gas_constant / temperature /
               temperature / temperature;
    }
    if (((variable1 == Variable::gas_phase_pressure) &&
         (variable2 == Variable::temperature)) ||
        ((variable1 == Variable::temperature) &&
         (variable2 == Variable::gas_phase_pressure)))
    {
        // d2rho_dpdT or d2rho_dTdp
        // extend to take pressure-temperature-dependent molar mass into account
        return -molar_mass / gas_constant / temperature / temperature;
    }

    OGS_FATAL(
        "IdealGasLaw::d2Value is implemented for derivatives with respect to "
        "phase pressure and temperature only.");

    return 0.;
}

}  // namespace MaterialPropertyLib
