/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CubicLawAfterShearSlip.h"
#include <cassert>

namespace
{
/// Stores the occurence of shear slip event.
struct CubicLawAfterShearSlipState
    : public MaterialLib::Fracture::Permeability::PermeabilityState
{
    bool shear_slip_occured = false;
};
}  // namespace

namespace MaterialLib::Fracture::Permeability
{
CubicLawAfterShearSlip::CubicLawAfterShearSlip(
    double const initial_creation_aperture,
    double const minimum_permeability,
    double const aperture_threshold)
    : _initial_creation_aperture(initial_creation_aperture),
      _minimum_permeability(minimum_permeability),
      _aperture_threshold(aperture_threshold)
{
}

void CubicLawAfterShearSlip::setShearSlipState(PermeabilityState& state,
                                               bool const value) const
{
    assert(dynamic_cast<CubicLawAfterShearSlipState*>(*state) != nullptr);
    static_cast<CubicLawAfterShearSlipState&>(state).shear_slip_occured = value;
}

bool CubicLawAfterShearSlip::shearSlipOccured(
    PermeabilityState const* const state) const
{
    assert(state != nullptr);
    return static_cast<CubicLawAfterShearSlipState const*>(state)
        ->shear_slip_occured;
}

double CubicLawAfterShearSlip::permeability(
    PermeabilityState const* const state,
    double const aperture0,
    double const aperture_m) const
{
    double const aperture = std::max(
        0.,
        shearSlipOccured(state)
            ? aperture_m + _initial_creation_aperture - _aperture_threshold
            : aperture0 - _aperture_threshold);

    return std::max(_minimum_permeability, aperture * aperture / 12);
}

double CubicLawAfterShearSlip::dpermeability_daperture(
    PermeabilityState const* const state,
    double const /*aperture0*/,
    double const aperture_m) const
{
    if (!shearSlipOccured(state))
    {
        return 0;
    }

    double const aperture =
        aperture_m + _initial_creation_aperture - _aperture_threshold;
    if (aperture <= 0)
    {
        return 0;
    }

    double const k = aperture * aperture / 12;
    if (k <= _minimum_permeability)
    {
        return 0;
    }

    return aperture / 6;
}

std::unique_ptr<PermeabilityState> CubicLawAfterShearSlip::getNewState() const
{
    return std::make_unique<CubicLawAfterShearSlipState>();
}
}  // namespace MaterialLib::Fracture::Permeability
