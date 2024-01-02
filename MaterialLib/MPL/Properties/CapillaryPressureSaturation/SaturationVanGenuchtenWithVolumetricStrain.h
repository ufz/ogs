/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
 * \brief A strain dependent bimodal water retention model.
 *
 * It is based on the van Genuchten model.
 *
 * The equation is as follows
 * \f[S_{e} = \left((e_{m}+a \Delta e)
 *   \left[\frac{1}{((\frac{p_{c}}{p_{b}})^{n}+1)}\right]^{m}
 *  + (e_{M} - a \Delta e)
 *   \left[\frac{1}{((\frac{p_{c}}{p_{{b,M}}})^{n} +1)}\right]^{n}
 * \right) \frac{1}{e}\f]
 * with effective saturation defined as \f$S_{e}=\frac{S-S_r}{S_{max}-S_r}\f$,
 * where \f$S_r\f$ and \f$S_{max}\f$ are the residual saturation and the maximum
 * saturation, respectively.
 * The exponent \f$m \in (0,1)\f$ is the same as in the van Genuchten equation (
 * see SaturationVanGenuchten). In the original work another exponent \f$n\f$ is
 * used, but usually set to \f$n = 1 / (1 - m)\f$, and also in this
 * implementation.
 * The pressure scaling parameter \f$p_{b}\f$ is added by the user and is the
 * scaling parameter of the micropores. The scaling parameter of the macropores
 * can be calculated as follows \f$p_{b,M} = p_{b} d_{diff}\f$ The total void
 * ratio and the void ratio of the micropores are \f$e_0\f$ and \f$e_m\f$.
 * Another scaling factor \f$a\f$ scales the effect of the strain
 *
 * The changing void ratio is calculated as
 * \f$\Delta e = -\frac{(1-e)\epsilon_{vol}}{e}\f$,
 * with \f$\epsilon_{vol}\f$ as the volumetric strain. The result is then
 * clamped between the residual and maximum liquid saturations.
 */
class SaturationVanGenuchtenWithVolumetricStrain final : public Property
{
public:
    SaturationVanGenuchtenWithVolumetricStrain(
        std::string name,
        double const residual_liquid_saturation,
        double const residual_gas_saturation,
        double const exponent,
        double const p_b,
        double const e_0,
        double const e_m,
        double const a,
        double const d_diff);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SaturationVanGenuchtenWithVolumetricStrain' is "
                "implemented on the 'media' scale only.");
        }
    }

    /// Those methods override the base class implementations and actually
    /// compute and set the property values_ and dValues_.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    double const S_L_res_;
    double const S_L_max_;
    double const m_;
    double const p_b_;
    double const e_0_;
    double const e_m_;
    double const a_;
    double const d_diff_;
};
}  // namespace MaterialPropertyLib
