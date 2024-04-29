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
struct DarcyVelocityData
{
    GlobalDimVector<DisplacementDim> w_GS;
    GlobalDimVector<DisplacementDim> w_LS;

    static auto reflect()
    {
        using Self = DarcyVelocityData<DisplacementDim>;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("velocity_gas", &Self::w_GS),
            R::makeReflectionData("velocity_liquid", &Self::w_LS)};
    }
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
