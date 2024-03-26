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
struct PorosityDerivativeData
{
    double dphi_S_dT = nan;
    double phi_0 = nan;  // Initial porosity.
};

struct PorosityData
{
    double phi = nan;

    static auto reflect()
    {
        using Self = PorosityData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{R::makeReflectionData("porosity", &Self::phi)};
    }
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
