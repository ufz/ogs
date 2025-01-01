/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
