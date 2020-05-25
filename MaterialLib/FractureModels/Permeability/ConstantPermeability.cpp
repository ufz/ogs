/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ConstantPermeability.h"

namespace MaterialLib::Fracture::Permeability
{
ConstantPermeability::ConstantPermeability(double const permeability)
    : permeability_(permeability)
{
}

double ConstantPermeability::permeability(
    PermeabilityState const* const /*state*/,
    double const /*aperture0*/,
    double const /*aperture_m*/) const
{
    return permeability_;
}

double ConstantPermeability::dpermeability_daperture(
    PermeabilityState const* const /*state*/,
    double const /*aperture0*/,
    double const /*aperture_m*/) const
{
    return 0.;
}
}  // namespace MaterialLib::Fracture::Permeability
