/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
 *   \brief The van Genuchten capillary pressure model:
 *
 *   \f[p_c(S)=p_b (S_e^{-1/m}-1)^{1-m}\f]
 *   with
 *   \f[S_e=\frac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &p_b&            \mbox{ entry pressure,}\\
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &m \in (0,1) &        \mbox{ exponent.}\\
 *    \f}
 *
 *    Note in some expressions, a parameter of \f$n\f$ is introduced, where
 *     \f[n=1/(1-m)\f].
 *
 *    If \f$\alpha\f$ instead of \f$p_b\f$ is available, \f$p_b\f$ can
 * be calculated
 * as
 *    \f[p_b=\rho g/\alpha\f].
 *
 *  If the capillary pressure is known, the saturation can be
 *  obtained by this model with
 *  \f[S(p_c)=\left \{
 *  \begin{array}{1}
 *   S_{\mbox{max}},\, p_c < 0,\\
 *    \left( \left(\dfrac{p_c}{p_b}\right)^{\frac{1}{1-m}} +1\right)^{-m}
 *     (S_{\mbox{max}}-S_r) +S_r,\, p_c  \geq 0
      \end{array}
 *   \right.
 * \f].
 *
 *   class SaturationVanGenuchten handles the computations associated
 *    with  \f$S(p_c)\f$,
 *  while class CapillaryPressureVanGenuchten and
 *  class CapillaryPressureRegularizedVanGenuchten deal with the
 *  computations of \f$ p_c(S) \f$.
 */
class SaturationVanGenuchten final : public Property
{
public:
    SaturationVanGenuchten(double const residual_liquid_saturation,
                           double const residual_gas_saturation,
                           double const exponent,
                           double const entry_pressure);

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

    /// gets \f$ S(p_c) \f$
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// gets \f$ \frac{\partial S(p_c)}{\partial  p_c} \f$
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

    /// gets \f$ \frac{\partial^2 S(p_c)}{\partial  p_c^2} \f$
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t,
                             double const dt) const override;

private:
    Medium* _medium = nullptr;
    double const _S_L_res;
    double const _S_L_max;
    double const _m;
    double const _p_b;
};
}  // namespace MaterialPropertyLib
