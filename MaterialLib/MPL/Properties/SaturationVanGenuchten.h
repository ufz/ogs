/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <cmath>
#include <limits>

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/**
 *   \brief van Genuchten water retention model
 *
 *   \f[p_c=p_b (S_e^{-1/m}-1)^{1-m}\f]
 *   with
 *   \f[S_e=\frac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &p_b&            \mbox{ entry pressure,}\\
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &m(<=1) &        \mbox{ exponent.}\\
 *    \f}
 *
 *    Note:
 *     \f[m=1/(1-n)\f].
 *
 *    If \f$\alpha\f$ instead of \f$p_b\f$ is available, \f$p_b\f$ can
 * be calculated
 * as
 *    \f[p_b=\rho g/\alpha\f]
 *
 * This property must be a medium property, it computes the saturation
 *  of the wetting phase as function of capillary pressure.
 */
class SaturationVanGenuchten final : public Property
{
private:
    Medium* _medium = nullptr;
    double const _S_L_res;
    double const _S_L_max;
    double const _m;
    double const _p_b;
    double const _pc_max; ///< Maximum capillary pressiure

public:
    SaturationVanGenuchten(double const residual_liquid_saturation,
                           double const residual_gas_saturation,
                           double const exponent,
                           double const entry_pressure,
                           double const max_capillary_pressure);

    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'SaturationVanGenuchten' is implemented on the "
                "'media' scale only.");
        }
        _medium = std::get<Medium*>(scale_pointer);
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property _values and _dValues.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;

    /// Gets \f$ \frac{\partial S}{\partial p_c} \f$
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;

    /// Gets \f$ \frac{\partial^2 S}{\partial p_c^2} \f$
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/,
                             double const /*dt*/) const override;

    /// This method gets capillary pressure via saturation.
    PropertyDataType inverse_value(VariableArray const& variable_array,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t,
                                   double const dt) const override;
};
}  // namespace MaterialPropertyLib
