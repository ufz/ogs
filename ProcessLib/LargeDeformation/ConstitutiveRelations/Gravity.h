// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/StrongType.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "SolidDensity.h"

namespace ProcessLib::LargeDeformation
{
template <int DisplacementDim>
using VolumetricBodyForce =
    BaseLib::StrongType<GlobalDimVector<DisplacementDim>, struct GravityTag>;

template <int DisplacementDim>
struct GravityModel
{
    explicit GravityModel(ProcessLib::ConstitutiveRelations::SpecificBodyForce<
                          DisplacementDim> const& specific_body_force)
        : specific_body_force_(specific_body_force)
    {
    }

    void eval(SolidDensity const& rho_SR,
              VolumetricBodyForce<DisplacementDim>& out) const;

    static GravityModel create(
        ProcessLib::ConstitutiveRelations::SpecificBodyForce<
            DisplacementDim> const& specific_body_force)
    {
        return GravityModel{specific_body_force};
    }

private:
    // TODO (naumov) Do we need to store this for each integration point?
    ProcessLib::ConstitutiveRelations::SpecificBodyForce<DisplacementDim> const
        specific_body_force_;
};

extern template struct GravityModel<2>;
extern template struct GravityModel<3>;
}  // namespace ProcessLib::LargeDeformation
