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
/// A constant permeability model.
class ConstantPermeability final : public Permeability
{
public:
    explicit ConstantPermeability(double const permeability);

private:
    double permeability(PermeabilityState const* const /*state*/,
                        double const /*aperture0*/,
                        double const /*aperture_m*/) const override;

    double dpermeability_daperture(PermeabilityState const* const /*state*/,
                                   double const /*aperture0*/,
                                   double const /*aperture_m*/) const override;

private:
    double const _permeability;
};
}  // namespace MaterialLib::Fracture::Permeability
