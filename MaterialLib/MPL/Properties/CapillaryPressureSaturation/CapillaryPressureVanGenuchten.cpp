/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 20, 2020, 9:59 AM
 */

#include "CapillaryPressureVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
void checkVanGenuchtenExponentRange(const double m)
{
    if (!(m > 0 && m < 1))
    {
        OGS_FATAL(
            "The exponent value m = %g of van Genuchten saturation model, is "
            "out of its range of (0, 1)",
            m);
    }
}

CapillaryPressureVanGenuchten::CapillaryPressureVanGenuchten(
    double const residual_liquid_saturation,
    double const residual_gas_saturation, double const exponent,
    double const entry_pressure, double const max_capillary_pressure)
    : _residual_saturation(residual_liquid_saturation),
      _maximuml_saturation(1. - residual_gas_saturation),
      _m(exponent),
      _p_b(entry_pressure),
      _pc_max(max_capillary_pressure)
{
    checkVanGenuchtenExponentRange(_m);
}

PropertyDataType CapillaryPressureVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S =
        std::clamp(saturation, _residual_saturation, _maximuml_saturation);
    const double Se = (S - _residual_saturation) /
                      (_maximuml_saturation - _residual_saturation);
    const double pc =
        (Se < 1.0) ? _p_b * std::pow(std::pow(Se, (-1.0 / _m)) - 1.0, 1.0 - _m)
                   : 0.0;
    return std::clamp(pc, 0.0, _pc_max);
}

PropertyDataType CapillaryPressureVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S =
        std::clamp(saturation, _residual_saturation, _maximuml_saturation);
    const double Se = (S - _residual_saturation) /
                      (_maximuml_saturation - _residual_saturation);

    if (!(Se < 1.0))
    {
        return 0.0;
    }
    const double val1 = std::pow(Se, -1.0 / _m);
    const double val2 = std::pow(val1 - 1.0, -_m);
    return _p_b * (_m - 1.0) * val1 * val2 / (_m * (S - _residual_saturation));
}

PropertyDataType CapillaryPressureVanGenuchten::d2Value(
    VariableArray const& variable_array, Variable const /*primary_variable1*/,
    Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S =
        std::clamp(saturation, _residual_saturation, _maximuml_saturation);
    const double Se = (S - _residual_saturation) /
                      (_maximuml_saturation - _residual_saturation);
    if (!(Se < 1.0))
    {
        return 0.0;
    }

    const double val1 = std::pow(Se, 1.0 / _m);
    return -_p_b /
           (_m * _m * (S - _residual_saturation) * (S - _residual_saturation)) *
           std::pow(1 - val1, -_m - 1) * std::pow(val1, _m - 1) *
           ((1 - _m * _m) * val1 + _m - 1);
}

CapillaryPressureRegularizedVanGenuchten::
    CapillaryPressureRegularizedVanGenuchten(
        double const residual_liquid_saturation,
        double const residual_gas_saturation,
        double const exponent,
        double const entry_pressure,
        double const max_capillary_pressure)
    : _saturation_nonwet_r(residual_gas_saturation),
      _residual_saturation(residual_liquid_saturation),
      _maximuml_saturation(1. - residual_gas_saturation),
      _m(exponent),
      _p_b(entry_pressure),
      _pc_max(max_capillary_pressure)
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

    double const Sg = 1 - saturation;
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

    double const Sg = 1 - saturation;
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
    VariableArray variable_array1;
    variable_array1[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] =
        saturation + perturbation;

    PropertyDataType val1 =
        dValue(variable_array1, primary_variable1, pos, t, dt);
    PropertyDataType val0 =
        dValue(variable_array, primary_variable1, pos, t, dt);

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
