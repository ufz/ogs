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
/// Permeability model including shear slip dependency. A minimum permeability
/// is maintained befor shear slip occurence. After a shear slip event the model
/// becomes a cubic law permeability model.
class CubicLawAfterShearSlip final : public Permeability
{
public:
    explicit CubicLawAfterShearSlip(double const initial_creation_aperture,
                                    double const minimum_permeability,
                                    double const aperture_threshold);

    void setShearSlipState(PermeabilityState& state, bool const value) const;

private:
    bool shearSlipOccured(PermeabilityState const* const state) const;

    double permeability(PermeabilityState const* const state,
                        double const aperture0,
                        double const aperture_m) const override;

    double dpermeability_daperture(PermeabilityState const* const state,
                                   double const /*aperture0*/,
                                   double const aperture_m) const override;

    std::unique_ptr<PermeabilityState> getNewState() const override;

private:
    /// An induced initial creation aperture. As soon as shear- or tensile-
    /// failure occur, an induced creation aperture is added.
    double const _initial_creation_aperture;

    double const _minimum_permeability;

    /// Aperture with no hydraulic effect.
    double const _aperture_threshold;
};
}  // namespace MaterialLib::Fracture::Permeability
