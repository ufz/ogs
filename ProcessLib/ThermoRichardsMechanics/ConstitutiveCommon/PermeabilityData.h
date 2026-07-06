// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct PermeabilityData
{
    double k_rel = 0;
    double dk_rel_dS_L = 0;
    GlobalDimMatrix<DisplacementDim> Ki;

    static auto reflect()
    {
        using Self = PermeabilityData<DisplacementDim>;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("intrinsic_permeability", &Self::Ki),
            R::makeReflectionData("relative_permeability", &Self::k_rel)};
    }
};
// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct PermeabilityData<2>;
extern template struct PermeabilityData<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
