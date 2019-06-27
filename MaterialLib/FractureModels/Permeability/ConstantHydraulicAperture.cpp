/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstantHydraulicAperture.h"

namespace MaterialLib::Fracture::Permeability
{
/// Hydraulic aperture equals the initial mechanical aperture s.t. permeability
/// is constant.
double ConstantHydraulicAperture::permeability(
    PermeabilityState const* const /*state*/,
    double const aperture0,
    double const /*aperture_m*/) const
{
    return aperture0 * aperture0 / 12;
}

double ConstantHydraulicAperture::dpermeability_daperture(
    PermeabilityState const* const /*state*/,
    double const /*aperture0*/,
    double const /*aperture_m*/) const
{
    return 0;
}
}  // namespace MaterialLib::Fracture::Permeability
