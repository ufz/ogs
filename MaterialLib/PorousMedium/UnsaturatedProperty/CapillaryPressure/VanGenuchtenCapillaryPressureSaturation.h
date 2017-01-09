/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   VanGenuchtenCapillaryPressureSaturation.h
 *
 *  Created on October 28, 2016, 6:05 PM
 */

#pragma once

#include <limits>

#include "CapillaryPressureSaturation.h"

namespace MaterialLib
{
namespace PorousMedium
{
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
 *    If \f$\alpha\f$ instead of \f$p_b\f$ is available, \f$p_b\f$ can be
 * calculated
 * as
 *    \f[p_b=\rho g/\alpha\f]
 */
class VanGenuchtenCapillaryPressureSaturation final
    : public CapillaryPressureSaturation
{
public:
    /**
     * @param pb     Entry pressure, \f$ p_b \f$
     * @param Sr     Residual saturation, \f$ S_r \f$
     * @param Smax   Maximum saturation, \f$ S_{\mbox{max}} \f$
     * @param m      Exponent, \f$ m \f$
     * @param Pc_max Maximum capillary pressure, \f$ P_c^{\mbox{max}}\f$
     */
    VanGenuchtenCapillaryPressureSaturation(const double pb, const double Sr,
                                            const double Smax, const double m,
                                            const double Pc_max)
        : CapillaryPressureSaturation(Sr, Smax, Pc_max), _pb(pb), _m(m)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "van Genuchten water retention model.";
    }

    /// Get capillary pressure.
    double getCapillaryPressure(const double saturation) const override;

    /// Get saturation.
    double getSaturation(const double capillary_pressure) const override;

    /// Get the derivative of the capillary pressure with respect to saturation
    double getdPcdS(const double saturation) const override;

private:
    const double _pb;  ///< Entry pressure.
    const double _m;   ///< Exponent m, m in [0,1]. n=1/(1-m).
};

}  // end namespace
}  // end namespace
