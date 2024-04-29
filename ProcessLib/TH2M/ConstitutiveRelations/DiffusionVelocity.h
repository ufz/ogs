/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct DiffusionVelocityData
{
    GlobalDimVector<DisplacementDim> d_CG;
    GlobalDimVector<DisplacementDim> d_WG;
    GlobalDimVector<DisplacementDim> d_CL;
    GlobalDimVector<DisplacementDim> d_WL;

    static auto reflect()
    {
        using Self = DiffusionVelocityData<DisplacementDim>;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("diffusion_velocity_gas_gas", &Self::d_CG),
            R::makeReflectionData("diffusion_velocity_vapour_gas", &Self::d_WG),
            R::makeReflectionData("diffusion_velocity_solute_liquid",
                                  &Self::d_CL),
            R::makeReflectionData("diffusion_velocity_liquid_liquid",
                                  &Self::d_WL)};
    }
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
