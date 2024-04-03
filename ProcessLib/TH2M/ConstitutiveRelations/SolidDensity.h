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
struct SolidDensityData
{
    double rho_SR = nan;

    static auto reflect()
    {
        using Self = SolidDensityData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("solid_density", &Self::rho_SR)};
    }
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
