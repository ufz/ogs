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
struct MassMoleFractionsData
{
    double xnCG = 0.;
    double xmCG = 0.;
    double xnWL = 0.;
    double xmWL = 0.;

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
