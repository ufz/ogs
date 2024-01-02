/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 1, 2020, 2:15 PM
 */

#pragma once

#include <limits>

#include "BaseLib/ConfigTree-fwd.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;

/**
 *  van Genuchten-Mualem relative permeability function for non-wetting phase
 *  in terms of effective wetting-phase saturation \cite lenhard1987model :
 *
 *   \f[k_{rel}^n= (1 - S_e)^{1/2} (1 - S_e^{1/m})^{2m}\f]
 *   with
 *   \f[S_e=\frac{S^L-S^L_r}{S^L_{\mbox{max}}-S^L_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &S^L_r&            \mbox{residual saturation of wetting phase,}\\
 *       &S^L_{\mbox{max}}& \mbox{maximum saturation of wetting phase,}\\
 *       &m\, \in (0, 1) &    \mbox{ exponent.}\\
 *    \f}
 *
 *  The derivative of the relative permeability with respect to saturation is
 *  computed as
 *  \f[\frac{\mathrm{d} k_{rel}^n}{\mathrm{d}S^L}=
 *  -(\dfrac{[1-S_e^{1/m}]^{2m}}{2\sqrt{1-S_e}}+2\sqrt{1-S_e}
 *  {(1-S_e^{1/m})}^{2*m-1}
 *   S_e^{1/m-1})/(S^L_\mbox{max}-S^L_r) \f]
 *   As \f$S^L \to S^L_\mbox{max}\f$, or \f$S_e \to 1\f$,
 *  \f$\dfrac{[1-S_e^{1/m}]^{2m}}{2\sqrt{1-S_e}}\f$ has a limit of zero.
 *
 */
class RelPermNonWettingPhaseVanGenuchtenMualem final : public Property
{
public:
    /**
     * @param name     Name of the property,
     * @param S_L_r    Residual saturation of the wetting phase,
     *                  \f$ S^L_r \f$
     * @param S_n_r    Residual saturation of the non-wetting phase,
     *                  \f$ S^n_{r} \f$
     * @param m        Exponent, \f$ m \in [0,1]\f$
     * @param krel_min Minimum relative permeability,
     *                  \f$ k_{rel}^n_{\mbox{min}}\f$
     */
    RelPermNonWettingPhaseVanGenuchtenMualem(std::string name,
                                             const double S_L_r,
                                             const double S_n_r,
                                             const double m,
                                             const double krel_min);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'RelPermNonWettingPhaseVanGenuchtenMualem' is "
                "implemented on the 'media' scale only.");
        }
    }

    /// \return \f$k_{rel}^n \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;

    /// \return \f$ \frac{\mathrm{d} k_{rel}^n}{\mathrm{d}S^L} \f$.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

    /**
     * Computes the saturation that gives the minimum relative permeability by
     * using the Regula–Falsi Method.
     * @return \f$ S^L\f$ that gives the minimum relative permeability.
     */
    double computeSaturationForMinimumRelativePermeability() const;

private:
    const double S_L_r_;     ///< Residual saturation of wetting phase.
    const double S_L_max_;   ///< Maximum saturation of wetting phase.
    const double m_;         ///< Exponent \f$ m \f$.
    const double krel_min_;  ///< Minimum relative permeability.
    const double
        S_L_for_krel_min_;  ///< Liquid saturation that gives \c krel_min_.
};
}  // namespace MaterialPropertyLib
