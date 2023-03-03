/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SaturationVolStrain.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include <limits>

#include "BaseLib/Error.h"
#include "MathLib/MathTools.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
SaturationVolStrain::SaturationVolStrain(
    std::string name,
    double const residual_liquid_saturation,
    double const residual_gas_saturation,
    double const exponent,
    double const p_b,
    double const b11,
    double const b22,
    double const b33)
    : S_L_res_(residual_liquid_saturation),
      S_L_max_(1. - residual_gas_saturation),
      m_(exponent),
      p_b_(p_b),
      b11_(b11),
      b22_(b22),
      b33_(b33)
{
    name_ = std::move(name);

    if (!(m_ > 0 && m_ < 1))
    {
        OGS_FATAL(
            "The exponent value m = {:g} of van Genuchten saturation model, is "
            "out of its range of (0, 1)",
            m_);
    }
}

PropertyDataType SaturationVolStrain::value(
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
    /*double const e_vol_pls = variable_array.equivalent_plastic_strain;*/

    double const n = 1. / (1. - m_);
    double const ten_base_exponent = e_vol > 0.0 ? b33_ * e_vol : b22_ * e_vol;
    double const factor = /*std::exp(b11_ * e_vol_pls)*/ std::pow(10.0, ten_base_exponent) ;

    double const p = (p_cap * factor) / p_b_;
    double const p_to_n = std::pow(p, n);
    double const S_eff = std::pow(p_to_n + 1., -m_);
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;

    return std::clamp(S, S_L_res_, S_L_max_);
}


PropertyDataType SaturationVolStrain::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{

    const double p_cap = variable_array.capillary_pressure;
    double const e_vol = variable_array.volumetric_strain;

    if (p_cap <= 0)
    {
        return 0.;
    }

    if (variable != Variable::capillary_pressure)
    {
        OGS_FATAL(
            "OrthotropicEmbeddedFracturePermeability::dValue is implemented "
            "for derivatives with respect to cap_p only.");
    }


    double const n = 1. / (1. - m_);
    double const ten_base_exponent = e_vol > 0.0 ? b33_ * e_vol : b22_ * e_vol;
    double const factor =
        std::pow(10.0, ten_base_exponent) /*std::exp(b11_ e_vol_pls)*/;
    double const factorn =
        std::pow(10.0, ten_base_exponent * n); 
    double const p_f = (p_cap * factor) / p_b_;
    double const p = p_cap / p_b_;
    double const p_to_n = std::pow(p, n);
    double const p_to_n_f = std::pow(p_f, n);
    double const S_eff = std::pow(p_to_n_f + 1., -m_);
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;

    if (S < S_L_res_ || S > S_L_max_)
    {
        return 0.;
    }

    switch (variable)
    {
        case Variable::capillary_pressure:
            {
                double const dS_eff_dp_cap = -m_ * n * factorn * p_to_n *
                                 std::pow(1 + factorn * p_to_n, -1. - m_) /
                                 (p_cap );
                return dS_eff_dp_cap * (S_L_max_ - S_L_res_);
            }
        case Variable::volumetric_strain:
            {
                double const dS_eff_de_vol = -m_ * n * std::log(10) * p_to_n * factorn * std::pow(p_to_n * factorn + 1, -m_ -1);
                return dS_eff_de_vol * (S_L_max_ - S_L_res_);
            }
        default:
            {
            OGS_FATAL(
            "WaterDensityIAPWSIF97Region1::dValue is implemented for derivatives with "
            "respect to capillary_pressure or volumetric_strain only.");
            }
    }

}

}  // namespace MaterialPropertyLib