/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 20, 2020, 9:59 AM
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
 * The the capillary pressure is computed from saturation as above but is cut
 * off at maximum capillary pressure given by user.
 */
class CapillaryPressureVanGenuchten : public Property
{
public:
    CapillaryPressureVanGenuchten(std::string name,
                                  double const residual_liquid_saturation,
                                  double const residual_gas_saturation,
                                  double const exponent,
                                  double const p_b,
                                  double const maximum_capillary_pressure);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'CapillaryPressureVanGenuchten' is implemented "
                "on the 'media' scale only.");
        }
    }

    /// \returns \f$ p_c(S) \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// \returns \f$ \frac{\partial p_c(S)}{\partial  S} \f$
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

private:
    double const S_L_res_;    ///< Residual saturation of liquid phase.
    double const S_L_max_;    ///< Maximum saturation of liquid phase.
    double const m_;          ///< Exponent.
    double const p_b_;        ///< Pressure scaling factor.
    double const p_cap_max_;  ///< Maximum capillary pressure.
};
}  // namespace MaterialPropertyLib
