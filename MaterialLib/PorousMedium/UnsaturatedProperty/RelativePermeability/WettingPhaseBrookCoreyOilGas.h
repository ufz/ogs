/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WettingPhaseBrookCoreyOilGas.h
 *
 * Created on November 1, 2016, 3:37 PM
 */

#ifndef OGS_WETTING_PHASE_BROOK_COREY_OIL_GAS_H
#define OGS_WETTING_PHASE_BROOK_COREY_OIL_GAS_H

#include <array>
#include <limits>  // std::numeric_limits

#include "RelativePermeability.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief BrookCorey oil-gas model: wetting phase
 *
 *   \f[k{rel}= S_e^{3 + 2 /m}\f]
 *   with
 *   \f[S_e=\dfrac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &m(>=1) &        \mbox{ exponent.}\\
 *    \f}
 */
class WettingPhaseBrookCoreyOilGas final : public RelativePermeability
{
public:
    /** \param parameters An array contains the five parameters:
     *                     [0] \f$ S_{r} \f$
     *                     [1] \f$ S_{\mbox{max}} \f$
     *                     [2] \f$ m \f$
     *                     [3] \f$ K_{rel}^{\mbox{min}}\f$
     */
    WettingPhaseBrookCoreyOilGas(std::array<double, 4> const& parameters)
        : _Sr(parameters[0]),
          _Smax(parameters[1]),
          _mm(parameters[2]),
          _Krel_min(parameters[3])
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Wetting phase Brook-Corey relative permeability model.";
    }

    /// Get relative permeability value.
    double getValue(const double saturation) const override;

private:
    const double _Sr;        ///< Residual saturation.
    const double _Smax;      ///< Maximum saturation.
    const double _mm;        ///< Exponent (>=1.0), n=1/(1-mm).
    const double _Krel_min;  ///< Minimum relative permeability
};

}  // end namespace
}  // end namespace
#endif /* OGS_WETTING_PHASE_BROOK_COREY_OIL_GAS_H */
