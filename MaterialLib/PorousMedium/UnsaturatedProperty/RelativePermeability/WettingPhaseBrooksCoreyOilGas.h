/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WettingPhaseBrooksCoreyOilGas.h
 *
 * Created on November 1, 2016, 3:37 PM
 */

#pragma once

#include "RelativePermeability.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief BrooksCorey oil-gas model: wetting phase
 *
 *   \f[k{rel}= S_e^{3 + 2 /m}\f]
 *   with
 *   \f[S_e=\frac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &m(>=1) &        \mbox{ exponent.}\\
 *    \f}
 */
class WettingPhaseBrooksCoreyOilGas final : public RelativePermeability
{
public:
    /**
     * @param Sr       Residual saturation, \f$ S_r \f$
     * @param Smax     Maximum saturation, \f$ S_{\mbox{max}} \f$
     * @param m        Exponent, \f$ m \f$
     * @param krel_min Minimum relative permeability,\f$ k_{rel}^{\mbox{min}}\f$
     */
    WettingPhaseBrooksCoreyOilGas(const double Sr, const double Smax,
                                  const double m, const double krel_min)
        : RelativePermeability(Sr, Smax), _m(m), _krel_min(krel_min)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Wetting phase Brooks-Corey relative permeability model.";
    }

    /// Get relative permeability value.
    double getValue(const double saturation) const override;

    /// Get the derivative of relative permeability with respect to saturation.
    /// \param saturation Wetting phase saturation
    double getdValue(const double saturation) const override;

private:
    const double _m;         ///< Exponent m, m>=1.0.
    const double _krel_min;  ///< Minimum relative permeability
};

}  // end namespace
}  // end namespace
