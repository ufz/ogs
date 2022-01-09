/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \brief The van Genuchten capillary pressure model.
 *
 * The van Genuchten capillary pressure model (\cite Genuchten1980) is:
 * \f[p_c(S)=p_b (S_\text{eff}^{-1/m}-1)^{1-m}\f]
 * with effective saturation defined as
 * \f[S_\text{eff}=\frac{S-S_r}{S_{\text{max}}-S_r}.\f]
 * Above, \f$S_r\f$ and \f$S_{\text{max}}\f$ are the residual and the maximum
 * saturations.
 * The exponent \f$m \in (0,1)\f$ and the pressure scaling parameter \f$p_b\f$
 * (it is equal to \f$\rho g/\alpha\f$ in original publication) are given by the
 * user.
 * The scaling parameter \f$p_b\f$ is given in same units as pressure.
 *
 * In the original work another exponent \f$n\f$ is used, but usually set to
 * \f$n = 1 / (1 - m)\f$, and also in this implementation.
 *
 * The saturation is computed from the capillary pressure as follows:
 * \f[S(p_c)=
 * \begin{cases}
 *  S_{\text{max}} & \text{for $p_c \leq 0$, and}\\
 *   \left( \left(\frac{p_c}{p_b}\right)^{\frac{1}{1-m}} +1\right)^{-m}
 *    (S_{\text{max}}-S_r) +S_r& \text{for $p_c > 0$.}
 * \end{cases}
 * \f]
 * The result is then clamped between the residual and maximum liquid
 * saturations.
 */
class SaturationVanGenuchten final : public Property
{
public:
    SaturationVanGenuchten(std::string name,
                           double const residual_liquid_saturation,
                           double const residual_gas_saturation,
                           double const exponent,
                           double const p_b);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SaturationVanGenuchten' is implemented on the "
                "'media' scale only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property values_ and dValues_.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/,
                             double const /*dt*/) const override;

private:
    double const S_L_res_;
    double const S_L_max_;
    double const m_;
    double const p_b_;
};
}  // namespace MaterialPropertyLib
