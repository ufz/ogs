/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearElasticIsotropic.h"

namespace MaterialLib
{
namespace Fracture
{
/// A penalty function for negative aperture supression used as a multiplier to
/// the normal fracture stiffness.
///
/// The derivative is continuous at aperture = aperture0 and aperture =
/// aperture_cutoff.
inline double logPenaltyDerivative(double const aperture0,
                                   double const aperture,
                                   double const aperture_cutoff)
{
    if (aperture >= aperture0)
        return 1;

    // Logarithmic penalty
    if (aperture > aperture_cutoff)
    {
        double const penalty = std::log(aperture / aperture0);
        return 1 + penalty * penalty +
               2 * penalty / aperture * (aperture - aperture0);
    }

    // Linear penalty below aperture cutoff
    {
        double const penalty = std::log(aperture_cutoff / aperture0);
        return 1 + penalty * penalty +
               2 * penalty / aperture_cutoff *
                   (2 * aperture - aperture_cutoff - aperture0);
    }
};

inline double logPenalty(double const aperture0,
                         double const aperture,
                         double const aperture_cutoff)
{
    if (aperture >= aperture0)
        return 1;

    // Logarithmic penalty
    if (aperture > aperture_cutoff)
    {
        double const penalty = std::log(aperture / aperture0);
        return 1 + penalty * penalty;
    }

    // Linear penalty below aperture cutoff
    {
        double const penalty = std::log(aperture_cutoff / aperture0);
        return 1 + penalty * penalty +
               2 * penalty / aperture_cutoff * (aperture - aperture_cutoff);
    }
};
}  // namespace Fracture
}  // namespace MaterialLib
