// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct MassMoleFractionsData
{
    double xnCG = nan;
    double xmCG = nan;
    double xnWL = nan;
    double xmWL = nan;

    static auto reflect()
    {
        using Self = MassMoleFractionsData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("mole_fraction_liquid", &Self::xnWL),
            R::makeReflectionData("mass_fraction_liquid", &Self::xmWL),
            R::makeReflectionData("mole_fraction_gas", &Self::xnCG),
            R::makeReflectionData("mass_fraction_gas", &Self::xmCG)};
    }
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
