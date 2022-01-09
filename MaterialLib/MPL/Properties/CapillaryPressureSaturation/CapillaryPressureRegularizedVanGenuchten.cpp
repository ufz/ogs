/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 20, 2020, 9:30 AM
 */

#include "CapillaryPressureRegularizedVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/CheckVanGenuchtenExponentRange.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
void checkSaturationRange(const double Sl)
{
    if (Sl < 0 || Sl > 1)
    {
        OGS_FATAL("The saturation of {:e} is out of its range of [0, 1]", Sl);
    }
}

CapillaryPressureRegularizedVanGenuchten::
    CapillaryPressureRegularizedVanGenuchten(
        double const residual_liquid_saturation,
        double const maximum_liquid_saturation,
        double const exponent,
        double const p_b)
    : Sg_r_(1.0 - maximum_liquid_saturation),
      Sg_max_(1.0 - residual_liquid_saturation),
      m_(exponent),
      p_b_(p_b),
      PcBarvGSg_Sg_max_(getPcBarvGSg(Sg_max_)),
      dPcdSvGBarSg_max_(getdPcdSvGBar(Sg_max_))
{
    checkVanGenuchtenExponentRange(m_);
}

PropertyDataType CapillaryPressureRegularizedVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double Sl = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    checkSaturationRange(Sl);

    double const Sg = 1 - Sl;
    if (!(Sg < Sg_r_ || Sg > Sg_max_))
    {
        return getPcBarvGSg(Sg);
    }
    if (Sg < Sg_r_)
    {
        return 0.0;
    }

    return PcBarvGSg_Sg_max_ + dPcdSvGBarSg_max_ * (Sg - Sg_max_);
}

PropertyDataType CapillaryPressureRegularizedVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (primary_variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "CapillaryPressureRegularizedVanGenuchten::dValue is implemented "
            "for derivatives with respect to liquid saturation only.");
    }

    const double Sl = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    checkSaturationRange(Sl);

    double const Sg = 1 - Sl;
    if (!(Sg < Sg_r_ || Sg > Sg_max_))
    {
        return -getdPcdSvGBar(Sg);
    }
    if (Sg < Sg_r_)
    {
        return 0.0;
    }

    return -dPcdSvGBarSg_max_;
}

double CapillaryPressureRegularizedVanGenuchten::getPcBarvGSg(
    double const Sg) const
{
    double const S_bar = getSBar(Sg);
    return getPcvGSg(S_bar) - getPcvGSg(Sg_r_ + (Sg_max_ - Sg_r_) * xi_ / 2);
}

double CapillaryPressureRegularizedVanGenuchten::getSBar(double const Sg) const
{
    return Sg_r_ + (1 - xi_) * (Sg - Sg_r_) + 0.5 * xi_ * (Sg_max_ - Sg_r_);
}

double CapillaryPressureRegularizedVanGenuchten::getPcvGSg(
    double const Sg) const
{
    double const Se = (Sg_max_ - Sg) / (Sg_max_ - Sg_r_);
    return p_b_ * std::pow(std::pow(Se, (-1.0 / m_)) - 1.0, 1.0 - m_);
}

double CapillaryPressureRegularizedVanGenuchten::getdPcdSvGBar(
    double const Sg) const
{
    double S_bar = getSBar(Sg);
    return getdPcdSvG(S_bar) * (1 - xi_);
}

double CapillaryPressureRegularizedVanGenuchten::getdPcdSvG(
    const double Sg) const
{
    double const n = 1 / (1 - m_);
    double const Se = (Sg_max_ - Sg) / (Sg_max_ - Sg_r_);
    auto const temp = std::pow(Se, (-1 / m_));
    return p_b_ * (1 / (m_ * n)) * (1 / (Sg_max_ - Sg_r_)) *
           std::pow(temp - 1, (1 / n) - 1) * temp / Se;
}

}  // namespace MaterialPropertyLib
