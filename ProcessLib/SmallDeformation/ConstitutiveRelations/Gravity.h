/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "BaseLib/StrongType.h"
#include "SolidDensity.h"

namespace ProcessLib::SmallDeformation
{
template <int DisplacementDim>
using VolumetricBodyForce =
    BaseLib::StrongType<GlobalDimVector<DisplacementDim>, struct GravityTag>;

template <int DisplacementDim>
struct GravityModel
{
    explicit GravityModel(
        Eigen::Vector<double, DisplacementDim> const& specific_body_force)
        : specific_body_force_(specific_body_force)
    {
    }

    void eval(SolidDensity const& rho_SR,
              VolumetricBodyForce<DisplacementDim>& out) const;

private:
    // TODO (naumov) Do we need to store this for each integration point?
    Eigen::Vector<double, DisplacementDim> const specific_body_force_;
};

extern template struct GravityModel<2>;
extern template struct GravityModel<3>;
}  // namespace ProcessLib::SmallDeformation
