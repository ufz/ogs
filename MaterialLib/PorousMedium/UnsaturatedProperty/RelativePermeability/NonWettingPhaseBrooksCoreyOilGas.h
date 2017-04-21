/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   NonWettingPhaseBrooksCoreyOilGas.h
 *
 * Created on November 2, 2016, 10:47 AM
 */

#pragma once

#include "RelativePermeability.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief BrooksCorey oil-gas model: non-wetting phase
 *
 *   \f[k{rel}= (1-S_e)^2  (1 - S_e^{1 + 2 /m})\f]
 *   with
 *   \f[S_e=\frac{S^w-S_r}{S^w_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &S^w_r&            \mbox{residual saturation of wetting phase,}\\
 *       &S^w_{\mbox{max}}& \mbox{maximum saturation of wetting phase,}\\
 *       &m(>=1) &        \mbox{ exponent.}\\
 *    \f}
 */
class NonWettingPhaseBrooksCoreyOilGas final : public RelativePermeability
{
public:
    /**
     * @param Snr       Residual saturation of the non-wetting phase,
     *                  \f$ S^n_r \f$
     * @param Snmax     Maximum saturation  of the non-wetting phase,
     *                  \f$ S^n_{\mbox{max}} \f$
     * @param m         Exponent, \f$ m \f$
     * @param krel_min  Minimum relative permeability,
     *                  \f$ k_{rel}^{\mbox{min}}\f$
     */
    NonWettingPhaseBrooksCoreyOilGas(const double Snr, const double Snmax,
                                     const double m, const double krel_min)
        : RelativePermeability(1. - Snmax, 1. - Snr), _m(m), _krel_min(krel_min)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Non-wetting phase Brook-Corey relative permeability model.";
    }

    /// Get relative permeability value.
    /// \param saturation_w Wetting phase saturation
    double getValue(const double saturation_w) const override;

    /// Get the derivative of relative permeability with respect to saturation.
    /// \param saturation_w Wetting phase saturation
    double getdValue(const double saturation_w) const override;

private:
    const double _m;         ///< Exponent m, m>=1.0.
    const double _krel_min;  ///< Minimum relative permeability
};

}  // end namespace
}  // end namespace
