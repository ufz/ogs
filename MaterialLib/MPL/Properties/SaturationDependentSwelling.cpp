/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/SaturationDependentSwelling.h"

#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/CoordinateSystem.h"

namespace MaterialPropertyLib
{
SaturationDependentSwelling::SaturationDependentSwelling(
    std::array<double, 3> swelling_pressures,
    std::array<double, 3>
        exponents,
    double const lower_saturation_limit,
    double const upper_saturation_limit,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : _p(std::move(swelling_pressures)),
      _lambda(std::move(exponents)),
      _S_min(lower_saturation_limit),
      _S_max(upper_saturation_limit),
      _local_coordinate_system(local_coordinate_system)
{
}

void SaturationDependentSwelling::setScale(
    std::variant<Medium*, Phase*, Component*> scale_pointer)
{
    if (std::holds_alternative<Phase*>(scale_pointer))
    {
        _phase = std::get<Phase*>(scale_pointer);
        if (_phase->name != "Solid")
        {
            OGS_FATAL(
                "The property 'SaturationDependentSwelling' must be "
                "given in the 'Solid' phase, not in '{:s}' phase.",
                _phase->name.c_str());
        }
    }
    else
    {
        OGS_FATAL(
            "The property 'SaturationDependentSwelling' is "
            "implemented on the 'phase' scales only.");
    }
}

PropertyDataType SaturationDependentSwelling::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const /*t*/,
    double const dt) const
{
    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    auto const S_L_dot = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation_rate)]);

    Eigen::Matrix<double, 3, 3> const e =
        _local_coordinate_system == nullptr
            ? Eigen::Matrix<double, 3, 3>::Identity()
            : _local_coordinate_system->transformation_3d(pos);

    Eigen::Matrix<double, 3, 3> delta_sigma_sw =
        Eigen::Matrix<double, 3, 3>::Zero();

    if (S_L < _S_min)
    {
        return delta_sigma_sw;   // still being zero.
    }

    double const S_L_prev = S_L - S_L_dot * dt;

    double const S_eff = std::clamp((S_L - _S_min) / (_S_max - _S_min), 0., 1.);
    double const S_eff_prev =
        std::clamp((S_L_prev - _S_min) / (_S_max - _S_min), 0., 1.);

    double const delta_S_eff = S_eff - S_eff_prev;

    // \Delta\sigma_{sw} = - \sum_i k_i (\lambda p S_{eff}^{(\lambda_i - 1)}
    // e_i \otimes e_i \Delta S_L / (S_{max} - S_{min}), where
    // e_i \otimes e_i is a square matrix with e_i,0^2 e_i,0*e_i,1 etc.
    for (int i = 0; i < 3; ++i)
    {
        Eigen::Matrix<double, 3, 3> const ei_otimes_ei =
            e.col(i) * e.col(i).transpose();

        delta_sigma_sw -=
            _lambda[i] * _p[i] * std::pow(S_eff, _lambda[i] - 1) * ei_otimes_ei;
    }

    return (delta_sigma_sw * delta_S_eff / dt).eval();
}

PropertyDataType SaturationDependentSwelling::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const /*t*/,
    double const dt) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_saturation) &&
           "SaturationDependentSwelling::dValue is implemented for "
           " derivatives with respect to liquid saturation only.");

    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    auto const S_L_dot = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation_rate)]);

    Eigen::Matrix<double, 3, 3> const e =
        _local_coordinate_system == nullptr
            ? Eigen::Matrix<double, 3, 3>::Identity()
            : _local_coordinate_system->transformation_3d(pos);

    Eigen::Matrix<double, 3, 3> delta_sigma_sw =
        Eigen::Matrix<double, 3, 3>::Zero();

    if (S_L < _S_min)
    {
        return delta_sigma_sw;   // still being zero.
    }

    double const S_L_prev = S_L - S_L_dot * dt;

    double const S_eff = std::clamp((S_L - _S_min) / (_S_max - _S_min), 0., 1.);
    double const S_eff_prev =
        std::clamp((S_L_prev - _S_min) / (_S_max - _S_min), 0., 1.);

    double const delta_S_eff = S_eff - S_eff_prev;

    // Heaviside(delta S_eff,sw)
    if (std::abs(delta_S_eff) <= 0)
    {
        return delta_sigma_sw;  // still being zero.
    }

    for (int i = 0; i < 3; ++i)
    {
        Eigen::Matrix<double, 3, 3> const ei_otimes_ei =
            e.col(i) * e.col(i).transpose();

        delta_sigma_sw +=
            _lambda[i] * _p[i] * std::pow(S_eff, _lambda[i] - 1) * ei_otimes_ei;
    }
    return (delta_sigma_sw / (_S_max - _S_min)).eval();
}
}  // namespace MaterialPropertyLib
