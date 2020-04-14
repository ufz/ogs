/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 1, 2020, 2:15 PM
 */

#pragma once

#include <limits>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 *  van Genuchten-Mualem relative permeability function for non-wetting phase
 *   \cite lenhard1987model
 *
 *   \f[k_{rel}^n= (1 - S_e)^{1/2} (1 - S_e^{1/m})^{2m}\f]
 *   with
 *   \f[S_e=\frac{S^w-S_r}{S^w_{\mbox{max}}-S^w_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &S^w_r&            \mbox{residual saturation of wetting phase,}\\
 *       &S^w_{\mbox{max}}& \mbox{maximum saturation of wetting phase,}\\
 *       &m\, \in (0, 1) &    \mbox{ exponent.}\\
 *    \f}
 */
class RelPermNonWettingPhaseVanGenuchtenMualem final : public Property
{
public:
    /**
    * @param Snr       Residual saturation of the non-wetting phase,
    *                  \f$ S^n_r \f$
    * @param Snmax     Maximum saturation  of the non-wetting phase,
    *                  \f$ S^n_{\mbox{max}} \f$
    * @param m         Exponent, \f$ m \in [0,1]\f$
    * @param krel_min  Minimum relative permeability,
    *                  \f$ k_{rel}^{\mbox{min}}\f$
    */
    RelPermNonWettingPhaseVanGenuchtenMualem(const double Snr,
                                             const double Snmax, const double m,
                                             const double krel_min);

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'RelPermNonWettingPhaseVanGenuchtenMualem' is "
                "implemented on the 'media' or 'Phase' scale only.");
        }
    }

    /// \return \f$k_{rel}^n \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;

    /// \return \f$ \dfrac{\partial k_{rel}^n}{\partial S^w} \f$.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    /// A small number for an offset to set the bound of S, the saturation, such
    /// that S in  [Sr+minor_offset_, Smax-minor_offset_].
    const double minor_offset_ = 1.e-9;

    const double S_L_res_;   ///< Residual saturation of wetting phase.
    const double S_L_max;    ///< Maximum saturation of wetting phase.
    const double m_;         ///< Exponent \f$ m \f$.
    const double krel_min_;  ///< Minimum relative permeability
};
}  // namespace MaterialPropertyLib
