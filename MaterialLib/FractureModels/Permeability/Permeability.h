/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>

namespace MaterialLib::Fracture::Permeability
{
struct PermeabilityState
{
    virtual ~PermeabilityState() = default;
};

/**
 * Interface for fracture permeability models.
 */
class Permeability
{
public:
    virtual double permeability(PermeabilityState const* const state,
                                double const aperture0,
                                double const aperture_m) const = 0;

    virtual double dpermeability_daperture(PermeabilityState const* const state,
                                           double const aperture0,
                                           double const aperture_m) const = 0;

    virtual ~Permeability() = default;

    virtual std::unique_ptr<PermeabilityState> getNewState() const
    {
        return nullptr;
    }
};
}  // namespace MaterialLib::Fracture::Permeability
