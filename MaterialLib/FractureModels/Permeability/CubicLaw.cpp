/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CubicLaw.h"

namespace MaterialLib::Fracture::Permeability
{
double CubicLaw::permeability(PermeabilityState const* const /*state*/,
                              double const /*aperture0*/,
                              double const aperture_m) const
{
    return aperture_m * aperture_m / 12;
}

double CubicLaw::dpermeability_daperture(
    PermeabilityState const* const /*state*/,
    double const /*aperture0*/,
    double const aperture_m) const
{
    return aperture_m / 6;
}
}  // namespace MaterialLib::Fracture::Permeability
