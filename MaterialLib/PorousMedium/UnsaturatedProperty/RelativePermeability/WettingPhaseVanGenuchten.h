/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WettingPhaseVanGenuchten.h
 *
 * Created on November 2, 2016, 11:24 AM
 */

#pragma once

#include "RelativePermeability.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief van Genuchten model model: wetting phase
 *
 *   \f[k{rel}=  \sqrt{S_e} (1-(1-S_e^{1/m})^m)^2)\f]
 *   with
 *   \f[S_e=\frac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &m(<=1) &        \mbox{ exponent.}\\
 *    \f}
 *
 *    Note:
 *     \f[m=1/(1-n)\f].
 */
class WettingPhaseVanGenuchten final : public RelativePermeability
{
public:
    /**
     * @param Sr       Residual saturation, \f$ S_r \f$
     * @param Smax     Maximum saturation, \f$ S_{\mbox{max}} \f$
     * @param m        Exponent, \f$ m \f$
     * @param krel_min Minimum relative permeability,\f$ k_{rel}^{\mbox{min}}\f$
     */
    WettingPhaseVanGenuchten(const double Sr, const double Smax, const double m,
                             const double krel_min)
        : RelativePermeability(Sr, Smax), _m(m), _krel_min(krel_min)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Wetting phase van Genuchten relative permeability model.";
    }

    /// Get relative permeability value.
    double getValue(const double saturation) const override;

    /// Get the derivative of relative permeability with respect to saturation.
    /// \param saturation Wetting phase saturation
    double getdValue(const double saturation) const override;

private:
    const double _m;         ///< Exponent m, m in [0, 1], n=1/(1-m).
    const double _krel_min;  ///< Minimum relative permeability
};

}  // end namespace
}  // end namespace
