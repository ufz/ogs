// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct PermeabilityData
{
    double k_rel;
    double dk_rel_dS_L;
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
}  // namespace ProcessLib::ThermoRichardsMechanics
