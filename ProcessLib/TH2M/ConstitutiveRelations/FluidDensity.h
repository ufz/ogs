// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct FluidDensityData
{
    // gas phase density
    double rho_GR = nan;

    // liquid phase density
    double rho_LR = nan;

    static auto reflect()
    {
        using Self = FluidDensityData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("gas_density", &Self::rho_GR),
            R::makeReflectionData("liquid_density", &Self::rho_LR)};
    }
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
