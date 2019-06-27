/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Permeability.h"

namespace MaterialLib::Fracture::Permeability
{
/// Hydraulic aperture equals the initial mechanical aperture s.t. permeability
/// is constant.
class ConstantHydraulicAperture final : public Permeability
{
    double permeability(PermeabilityState const* const /*state*/,
                        double const aperture0,
                        double const /*aperture_m*/) const override;

    double dpermeability_daperture(PermeabilityState const* const /*state*/,
                                   double const /*aperture0*/,
                                   double const /*aperture_m*/) const override;
};
}  // namespace MaterialLib::Fracture::Permeability
