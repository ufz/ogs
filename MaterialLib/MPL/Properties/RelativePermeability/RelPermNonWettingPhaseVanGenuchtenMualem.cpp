/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 1, 2020, 2:15 PM
 */

#include "RelPermNonWettingPhaseVanGenuchtenMualem.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/CheckVanGenuchtenExponentRange.h"
#include "MathLib/Nonlinear/Root1D.h"

namespace MaterialPropertyLib
{
double computeVanGenuchtenMualemValue(const double S_L, const double S_L_r,
                                      const double S_L_max, const double m)
{
    const double Se = (S_L - S_L_r) / (S_L_max - S_L_r);
    return std::sqrt(1.0 - Se) * std::pow(1.0 - std::pow(Se, 1.0 / m), 2.0 * m);
}

RelPermNonWettingPhaseVanGenuchtenMualem::
    RelPermNonWettingPhaseVanGenuchtenMualem(std::string name,
                                             const double S_L_r,
                                             const double S_n_r,
                                             const double m,
                                             const double krel_min)
    : S_L_r_(S_L_r),
      S_L_max_(1. - S_n_r),
      m_(m),
      krel_min_(krel_min),
      S_L_for_krel_min_(computeSaturationForMinimumRelativePermeability())
{
    name_ = std::move(name);
    checkVanGenuchtenExponentRange(m_);
}

PropertyDataType RelPermNonWettingPhaseVanGenuchtenMualem::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_L = std::clamp(
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation)]),
        S_L_r_, S_L_max_);

    const double krel =
        computeVanGenuchtenMualemValue(S_L, S_L_r_, S_L_max_, m_);
    return std::max(krel_min_, krel);
}

PropertyDataType RelPermNonWettingPhaseVanGenuchtenMualem::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "RelPermNonWettingPhaseVanGenuchtenMualem::dValue is implemented "
            "for the derivative with respect to liquid saturation only.");
    }

    const double S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    if (S_L < S_L_r_ || S_L > S_L_for_krel_min_)
    {
        return 0.0;
    }

    if (std::fabs(S_L - S_L_max_) < std::numeric_limits<double>::epsilon())
    {
        return 0.0;
    }

    const double Se = (S_L - S_L_r_) / (S_L_max_ - S_L_r_);

    const double val1 = std::sqrt(1.0 - Se);
    const double val2 = 1.0 - std::pow(Se, 1.0 / m_);

    return (-0.5 * std::pow(val2, 2.0 * m_) / val1 -
            2.0 * std::pow(Se, -1.0 + 1.0 / m_) * val1 *
                std::pow(val2, 2.0 * m_ - 1.0)) /
           (S_L_max_ - S_L_r_);
}

double RelPermNonWettingPhaseVanGenuchtenMualem::
    computeSaturationForMinimumRelativePermeability() const
{
    auto f = [this](double S_L) -> double
    {
        return computeVanGenuchtenMualemValue(S_L, this->S_L_r_, this->S_L_max_,
                                              this->m_) -
               this->krel_min_;
    };

    auto root_finder =
        MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(
            f, S_L_r_, S_L_max_);

    root_finder.step(1000);

    return root_finder.getResult();
}
}  // namespace MaterialPropertyLib
