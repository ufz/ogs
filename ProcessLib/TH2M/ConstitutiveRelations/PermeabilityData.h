// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    double k_rel_G = nan;
    double k_rel_L = nan;
    double dk_rel_G_dS_L = nan;
    double dk_rel_L_dS_L = nan;
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
