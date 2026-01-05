// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PengRobinson.h"

#include <algorithm>
#include <boost/math/special_functions/pow.hpp>
#include <cmath>

#include "MathLib/Nonlinear/CubicRoots.h"

namespace MaterialPropertyLib
{

PropertyDataType PengRobinson::value(
    MaterialPropertyLib::VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double gas_constant = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pressure = variable_array.gas_phase_pressure;
    const double temperature = variable_array.temperature;
    const double molar_mass = variable_array.molar_mass;

    const double Tr = temperature / Tc_;  // reduced Temperature

    const double kappa = 0.37464 + 1.5422 * omega_ - 0.26992 * omega_ * omega_;

    const double alpha = boost::math::pow<2>(1 + kappa * (1 - std::sqrt(Tr)));

    // EOS in the form: 0 = rho^3 + z1*rho^2 + z2*rho + z3
    const double denominator =
        b_ *
        (pressure * b_ * b_ + b_ * gas_constant * temperature - a_ * alpha);

    const double z1 =
        (molar_mass * a_ * alpha - 3 * molar_mass * b_ * b_ * pressure -
         2 * molar_mass * gas_constant * temperature * b_) /
        denominator;
    const auto z2 = (molar_mass * molar_mass *
                     (b_ * pressure - gas_constant * temperature)) /
                    denominator;
    const auto z3 =
        (molar_mass * molar_mass * molar_mass * pressure) / denominator;

    MathLib::CubicSolver cubic_solver_(1., z1, z2, z3);
    return cubic_solver_.smallestPositiveRealRoot();
}

PropertyDataType PengRobinson::dValue(
    MaterialPropertyLib::VariableArray const& variable_array,
    Variable const variable, ParameterLib::SpatialPosition const& pos,
    double const t, double const dt) const
{
    if (variable == Variable::temperature)
    {
        const double temperature = variable_array.temperature;
        const double epsilon = 1.e-6;

        MaterialPropertyLib::VariableArray perturbed = variable_array;
        // Increase temperature by +epsilon
        perturbed.temperature = temperature + epsilon;
        const double rho_plus = std::get<double>(value(perturbed, pos, t, dt));

        // Decrease temperature by -epsilon
        perturbed.temperature = temperature - epsilon;
        const double rho_minus = std::get<double>(value(perturbed, pos, t, dt));

        // Calculate the central difference quotient
        return (rho_plus - rho_minus) / (2 * epsilon);
    }

    if (variable == Variable::gas_phase_pressure)
    {
        const double pressure = variable_array.gas_phase_pressure;
        const double epsilon = 1.e-3;

        MaterialPropertyLib::VariableArray perturbed = variable_array;
        // Increase pressure by +epsilon
        perturbed.gas_phase_pressure = pressure + epsilon;
        const double rho_plus = std::get<double>(value(perturbed, pos, t, dt));

        // Decrease pressure by -epsilon
        perturbed.gas_phase_pressure = pressure - epsilon;
        const double rho_minus = std::get<double>(value(perturbed, pos, t, dt));

        // Calculate the central difference quotient
        return (rho_plus - rho_minus) / (2 * epsilon);
    }

    OGS_FATAL(
        "PengRobinson::dValue is implemented for derivatives with respect to "
        "gas phase pressure or temperature only.");

    return 0.;
}

}  // namespace MaterialPropertyLib
