/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   NonWettingPhaseVanGenuchten.h
 *
 * Created on November 2, 2016, 11:24 AM
 */
#ifndef OGS_NON_WETTING_PHASE_VAN_GENUCHTEN_H
#define OGS_NON_WETTING_PHASE_VAN_GENUCHTEN_H

#include <array>
#include <limits>  // std::numeric_limits

#include "RelativePermeability.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief BrookCorey oil-gas model: non-wetting phase
 *
 *   \f[k{rel}= (1 - S_e)^{1/3} (1 - S_e^{1/m})^{2m}\f]
 *   with
 *   \f[S_e=\dfrac{S-S_r}{S_{\mbox{max}}-S_r}\f]
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
class NonWettingPhaseVanGenuchten final : public RelativePermeability
{
public:
    /** \param parameters An array contains the five parameters:
     *                     [0] \f$ S_{nr} \f$
     *                     [1] \f$ S_{n\mbox{max}} \f$
     *                     [2] \f$ m \f$
     *                     [3] \f$ K_{rel}^{\mbox{min}}\f$
     */
    NonWettingPhaseVanGenuchten(std::array<double, 4> const& parameters)
        : _Sr(1. - parameters[1]),
          _Smax(1. - parameters[0]),
          _mm(parameters[2]),
          _Krel_min(parameters[3])
    {
    }

    NonWettingPhaseVanGenuchten(const NonWettingPhaseVanGenuchten& orig)
        : _Sr(orig._Sr),
          _Smax(orig._Smax),
          _mm(orig._mm),
          _Krel_min(orig._Krel_min)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Non-wetting phase van Genuchten model.";
    }

    /// Get relative permeability value.
    /// \param saturation_w Non-wetting phase saturation
    double getValue(const double saturation_w) const override;

private:
    const double _Sr;        ///< Residual saturation of wetting phase, 1-Snr.
    const double _Smax;      ///< Maximum saturation of wetting phase., 1-Sn_max
    const double _mm;        ///< Exponent (<=1.0), n=1/(1-mm).
    const double _Krel_min;  ///< Minimum relative permeability

    const double _perturbation = std::numeric_limits<double>::epsilon();
};

}  // end namespace
}  // end namespace
#endif /* OGS_NON_WETTING_PHASE_VAN_GENUCHTEN_H */
