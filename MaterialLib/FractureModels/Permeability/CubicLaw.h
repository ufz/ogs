/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Permeability.h"

namespace MaterialLib::Fracture::Permeability
{
/// Hydraulic aperture equals the mechanical aperture s.t. multiplication of the
/// permeability by the mechanical aperture yields the cubic law.
class CubicLaw final : public Permeability
{
    double permeability(PermeabilityState const* const /*state*/,
                        double const /*aperture0*/,
                        double const aperture_m) const override;

    double dpermeability_daperture(PermeabilityState const* const /*state*/,
                                   double const /*aperture0*/,
                                   double const aperture_m) const override;
};
}  // namespace MaterialLib::Fracture::Permeability
