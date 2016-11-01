/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   BrookCoreyCapillaryPressure.h
 *
 *  Created on November 1, 2016, 9:45 AM
 */

#ifndef OGS_BROOK_COREY_CAPILLARY_PRESSURE_H
#define OGS_BROOK_COREY_CAPILLARY_PRESSURE_H

#include <array>
#include <limits>  // std::numeric_limits

#include "CapillaryPressureSaturation.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief Brook-Corey capillary pressure saturation model
 *
 *   \f[p_c=p_b S_e^{-1/m}\f]
 *   with
 *   \f[S_e=\dfrac{S-S_r}{S_{\mbox{max}}-S_r}\f]
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
class BrookCoreyCapillaryPressure final : public CapillaryPressureSaturation
{
public:
    /** \param parameters An array contains the five parameters:
     *                     [0] $f p_b $f
     *                     [1] $f S_r $f
     *                     [2] $f S_{\mbox{max}} $f
     *                     [3] $f m $f
     *                     [4] $f P_c^{\mbox{max}}$f
     */
    BrookCoreyCapillaryPressure(std::array<double, 5> const& parameters)
        : _pb(parameters[0]),
          _Sr(parameters[1]),
          _Smax(parameters[2]),
          _mm(parameters[3]),
          _Pc_max(parameters[4])
    {
    }

    BrookCoreyCapillaryPressure(const BrookCoreyCapillaryPressure& orig)
        : _pb(orig._pb),
          _Sr(orig._Sr),
          _Smax(orig._Smax),
          _mm(orig._mm),
          _Pc_max(orig._Pc_max)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Brook-Corey capillary pressure saturation model.";
    }

    /// Get capillary pressure.
    double getCapillaryPressure(const double saturation) const override;

    /// Get capillary pressure.
    double getSturation(const double capillary_pressure) const override;

    /// Get the derivative of the capillary pressure with respect to saturation
    double getdPcdS(const double saturation) const override;

private:
    const double _pb;      ///< Entry pressure.
    const double _Sr;      ///< Residual saturation.
    const double _Smax;    ///< Maximum saturation.
    const double _mm;      ///< Exponent (<=1.0), n=1/(1-mm).
    const double _Pc_max;  ///< Maximum capillaray pressure

    const double _perturbation = std::numeric_limits<double>::epsilon();
};

}  // end namespace
}  // end namespace
#endif /* OGS_BROOK_COREY_CAPILLARY_PRESSURE_H */
