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
struct PermeabilityData
{
    double k_rel_G;
    double k_rel_L;
    double dk_rel_G_dS_L;
    double dk_rel_L_dS_L;
    GlobalDimMatrix<DisplacementDim> Ki;

    static auto reflect()
    {
        using Self = PermeabilityData<DisplacementDim>;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("intrinsic_permeability", &Self::Ki),
            R::makeReflectionData("relative_permeability_gas", &Self::k_rel_G),
            R::makeReflectionData("relative_permeability_liquid",
                                  &Self::k_rel_L)};
    }
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
