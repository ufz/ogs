/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SaturationVanGenuchtenWithVolumetricStrain.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/MathTools.h"

namespace MaterialPropertyLib
{
SaturationVanGenuchtenWithVolumetricStrain::
    SaturationVanGenuchtenWithVolumetricStrain(
        std::string name,
        double const residual_liquid_saturation,
        double const residual_gas_saturation,
        double const exponent,
        double const p_b,
        double const e_0,
        double const e_m,
        double const a,
        double const d_diff)
    : S_L_res_(residual_liquid_saturation),
      S_L_max_(1. - residual_gas_saturation),
      m_(exponent),
      p_b_(p_b),
      e_0_(e_0),
      e_m_(e_m),
      a_(a),
      d_diff_(d_diff)
{
    name_ = std::move(name);

    if (!(m_ > 0 && m_ < 1))
    {
        OGS_FATAL("The exponent value m = {:g}, is out of its range of (0, 1)",
                  m_);
    }
}

PropertyDataType SaturationVanGenuchtenWithVolumetricStrain::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double p_cap = variable_array.capillary_pressure;

    if (p_cap <= 0)
    {
        return S_L_max_;
    }

    double const e_vol = variable_array.volumetric_strain;

    double const n = 1. / (1. - m_);
    double const d_e = -1 * (1 + e_0_) * e_vol / e_0_;
    double const p_b_M = p_b_ * (1 / d_diff_);
    double const p = p_cap / p_b_;
    double const p_to_n = std::pow(p, n);
    double const S_eff_mi = std::pow(p_to_n + 1., -m_);
    double const p_M = p_cap / p_b_M;
    double const p_to_n_M = std::pow(p_M, n);
    double const S_eff_M = std::pow(p_to_n_M + 1., -m_);
    double const S_eff =
        (e_m_ * S_eff_mi + ((e_0_ - e_m_ - (a_ * d_e)) * S_eff_M)) /
        (e_0_ - (a_ * d_e));
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;

    return std::clamp(S, S_L_res_, S_L_max_);
}

PropertyDataType SaturationVanGenuchtenWithVolumetricStrain::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::capillary_pressure)
    {
        OGS_FATAL(
            "SaturationVanGenuchtenWithVolumetricStrain::dValue is implemented "
            "for derivatives with respect to capillary pressure only.");
    }
    const double p_cap = variable_array.capillary_pressure;

    if (p_cap <= 0)
    {
        return 0.;
    }

    double const e_vol = variable_array.volumetric_strain;

    double const n = 1. / (1. - m_);
    double const d_e = -1 * (1 + e_0_) * e_vol / e_0_;
    double const p_b_M = p_b_ * (1 / d_diff_);
    double const p = p_cap / p_b_;
    double const p_to_n = std::pow(p, n);
    double const S_eff_mi = std::pow(p_to_n + 1., -m_);
    double const p_M = p_cap / p_b_M;
    double const p_to_n_M = std::pow(p_M, n);
    double const S_eff_M = std::pow(p_to_n_M + 1., -m_);
    double const S_eff =
        (e_m_ * S_eff_mi + ((e_0_ - e_m_ - (a_ * d_e)) * S_eff_M)) /
        (e_0_ - (a_ * d_e));
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;

    if (S < S_L_res_ || S > S_L_max_)
    {
        return 0.;
    }

    double const dS_eff_dp_cap =
        (-(e_m_)*n * m_ * p_to_n * std::pow(p_to_n + 1., -m_ - 1) -
         (e_0_ - e_m_ - (a_ * d_e)) * n * m_ * p_to_n_M *
             std::pow(p_to_n_M + 1., -m_ - 1)) /
        ((e_0_ - (a_ * d_e)) * p_cap);
    return dS_eff_dp_cap * (S_L_max_ - S_L_res_);
}
}  // namespace MaterialPropertyLib
