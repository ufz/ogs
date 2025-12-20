// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct VapourPartialPressureData
{
    // water partial pressure in gas phase
    double pWGR = nan;

    static auto reflect()
    {
        using Self = VapourPartialPressureData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("vapour_pressure", &Self::pWGR)};
    }
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
