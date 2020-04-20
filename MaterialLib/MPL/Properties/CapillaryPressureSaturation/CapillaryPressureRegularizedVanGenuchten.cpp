/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 20, 2020, 9:30 AM
 */

#include "CapillaryPressureRegularizedVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/CheckVanGenuchtenExponentRange.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
CapillaryPressureRegularizedVanGenuchten::
    CapillaryPressureRegularizedVanGenuchten(
        double const residual_liquid_saturation,
        double const maximum_liquid_saturation,
        double const exponent,
        double const p_b)
    : _saturation_nonwet_r(1.0- maximum_liquid_saturation),
      _residual_saturation(residual_liquid_saturation),
      _maximuml_saturation(maximum_liquid_saturation),
      _m(exponent),
      _p_b(p_b)
{
    checkVanGenuchtenExponentRange(_m);
}

PropertyDataType CapillaryPressureRegularizedVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    const double Sl =
        std::clamp(saturation, _residual_saturation, _maximuml_saturation);

    double const Sg = 1 - Sl;
    if (!(Sg < _saturation_nonwet_r || Sg > 1 - _residual_saturation))
    {
        return getPcBarvGSg(Sg);
    }
    if (Sg < _saturation_nonwet_r)
    {
        return getPcBarvGSg(_saturation_nonwet_r) +
               getdPcdSvGBar(_saturation_nonwet_r) *
                   (Sg - _saturation_nonwet_r);
    }

    return getPcBarvGSg(1 - _residual_saturation) +
           getdPcdSvGBar(1 - _residual_saturation) *
               (Sg - 1 + _residual_saturation);
}

PropertyDataType CapillaryPressureRegularizedVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    const double Sl =
        std::clamp(saturation, _residual_saturation, _maximuml_saturation);

    double const Sg = 1 - Sl;
    if (Sg >= _saturation_nonwet_r && Sg <= 1 - _residual_saturation)
    {
        return -getdPcdSvGBar(Sg);
    }
    if (Sg < _saturation_nonwet_r)
    {
        return -getdPcdSvGBar(_saturation_nonwet_r);
    }

    return -getdPcdSvGBar(1 - _residual_saturation);
}

PropertyDataType CapillaryPressureRegularizedVanGenuchten::d2Value(
    VariableArray const& variable_array, Variable const primary_variable1,
    Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double perturbation =
        std::sqrt(std::numeric_limits<double>::epsilon());

    VariableArray variable_array0;
    const double S0 = (std::fabs(saturation - _maximuml_saturation) <
                       std::numeric_limits<double>::epsilon())
                          ? saturation - perturbation
                          : saturation;
    variable_array0[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = S0;

    VariableArray variable_array1;
    variable_array1[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = S0 + perturbation;

    PropertyDataType val1 =
        dValue(variable_array1, primary_variable1, pos, t, dt);
    PropertyDataType val0 =
        dValue(variable_array0, primary_variable1, pos, t, dt);

    return (std::get<double>(val1) - std::get<double>(val0)) / perturbation;
}

double CapillaryPressureRegularizedVanGenuchten::getPcBarvGSg(
    double const Sg) const
{
    double const Sg_r = _saturation_nonwet_r;
    double const S_lr = _residual_saturation;
    double const S_bar = getSBar(Sg);
    return getPcvGSg(S_bar) - getPcvGSg(Sg_r + (1 - Sg_r - S_lr) * _xi / 2);
}

double CapillaryPressureRegularizedVanGenuchten::getSBar(double const Sg) const
{
    double const Sg_r = _saturation_nonwet_r;
    double const S_lr = _residual_saturation;
    return Sg_r + (1 - _xi) * (Sg - Sg_r) + 0.5 * _xi * (1 - Sg_r - S_lr);
}

double CapillaryPressureRegularizedVanGenuchten::getPcvGSg(
    double const Sg) const
{
    double const Sg_r = _saturation_nonwet_r;
    double const S_lr = _residual_saturation;
    double const S_le = (1 - Sg - S_lr) / (1 - Sg_r - _residual_saturation);
    return _p_b * std::pow(std::pow(S_le, (-1.0 / _m)) - 1.0, 1.0 - _m);
}

double CapillaryPressureRegularizedVanGenuchten::getdPcdSvGBar(
    double const Sg) const
{
    double S_bar = getSBar(Sg);
    return getdPcdSvG(S_bar) * (1 - _xi);
}

double CapillaryPressureRegularizedVanGenuchten::getdPcdSvG(
    const double Sg) const
{
    double const Sg_r = _saturation_nonwet_r;
    double const S_lr = _residual_saturation;
    double const n = 1 / (1 - _m);
    double const S_le = (1 - Sg - S_lr) / (1 - Sg_r - S_lr);
    return _p_b * (1 / (_m * n)) * (1 / (1 - S_lr - Sg_r)) *
           std::pow(std::pow(S_le, (-1 / _m)) - 1, (1 / n) - 1) *
           std::pow(S_le, (-1 / _m)) / S_le;
}

}  // namespace MaterialPropertyLib
