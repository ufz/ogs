/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   BrooksCoreyCapillaryPressureSaturation.h
 *
 *  Created on November 1, 2016, 9:45 AM
 */

#pragma once

#include <limits>  // std::numeric_limits

#include "CapillaryPressureSaturation.h"

#include "BaseLib/Error.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief Brooks-Corey capillary pressure saturation model
 *
 *   \f[p_c=p_b S_e^{-1/m}\f]
 *   with
 *   \f[S_e=\frac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &p_b&            \mbox{ entry pressure,}\\
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &m(>=1) &        \mbox{ exponent.}\\
 *    \f}
 *
 *    Note:
 *     \f[m=1/(1-n)\f]
 */
class BrooksCoreyCapillaryPressureSaturation final
    : public CapillaryPressureSaturation
{
public:
    /**
     * @param pb     Entry pressure, \f$ p_b \f$
     * @param Sr     Residual saturation, \f$ S_r \f$
     * @param Sg_r     Residual saturation of gas phase, \f$ Sg_r \f$
     * @param Smax   Maximum saturation, \f$ S_{\mbox{max}} \f$
     * @param m      Exponent, \f$ m \f$
     * @param Pc_max Maximum capillary pressure, \f$ P_c^{\mbox{max}}\f$
     */
    BrooksCoreyCapillaryPressureSaturation(const double pb, const double Sr,
                                           const double Sg_r, const double Smax,
                                           const double m, const double Pc_max)
        : CapillaryPressureSaturation(Sr, Sg_r, Smax, Pc_max), _pb(pb), _m(m)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Brooks-Corey capillary pressure saturation model.";
    }

    /// Get capillary pressure.
    double getCapillaryPressure(const double saturation) const override;

    /// Get saturation.
    double getSaturation(const double capillary_pressure) const override;

    /// Get the derivative of the capillary pressure with respect to saturation
    double getdPcdS(const double saturation) const override;

    /// Get the second derivative of the capillary pressure with respect to
    /// saturation
    double getd2PcdS2(const double /*saturation*/) const override
    {
        OGS_FATAL(
            "Second derivative of BrooksCorey d^2p_cap / dS^2 is not "
            "implemented.");
    }

private:
    const double _pb;  ///< Entry pressure.
    const double _m;   ///< Exponent m, m>1.
};

}  // end namespace
}  // end namespace
