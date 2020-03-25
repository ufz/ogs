/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <limits>
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 *   \brief The Brooks-Corey capillary pressure saturation model:
 *
 *   \f[p_c=p_b S_e^{-1/\lambda}\f]
 *   with
 *   \f[S_e=\frac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &p_b&            \mbox{ entry pressure,}\\
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &\lambda \in [0,1) &        \mbox{ exponent.}\\
 *    \f}
 *
 *  If the capillary pressure is known, the saturation can be
 *  obtained by this model with
 *  \f[
 *   S(p_c)=
 *   \begin{cases}
 *     S_{\mbox{max}}, & p_c < p_b,\\
 *    \left(\dfrac{p_c}{p_b}\right)^{-\lambda}
 *     (S_{\mbox{max}}-S_r) +S_r,& p_c  \geq p_b
 *  \end{cases}
 * \f].
 *
 *   class SaturationBrooksCorey handles the computations associated
 *    with  \f$S(p_c)\f$,
 *  while class CapillaryPressureVanGenuchten deals with the
 *  computations of
 *  \f$ p_c(S) \f$.
 *
 */
class SaturationBrooksCorey : public Property
{
public:
    SaturationBrooksCorey(const double residual_liquid_saturation,
                          const double residual_gas_saturation,
                          const double exponent,
                          const double entry_pressure);

    void setScale(std::variant<Medium*, Phase*, Component*> scale_pointer)
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'SaturationBrooksCorey' is implemented on the "
                "'media' scale only.");
        }
        _medium = std::get<Medium*>(scale_pointer);
    }

    /// gets \f$ S(p_c) \f$
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const;

    /// gets \f$ \frac{\partial S(p_c)}{\partial  p_c} \f$
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const;

    /// gets \f$ \frac{\partial^2 S(p_c)}{\partial  p_c^2} \f$
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t, double const dt) const;

protected:
    Medium* _medium = nullptr;
    const double _residual_liquid_saturation;
    const double _max_liquid_saturation;
    const double _exponent;
    const double _entry_pressure;
};
}  // namespace MaterialPropertyLib
