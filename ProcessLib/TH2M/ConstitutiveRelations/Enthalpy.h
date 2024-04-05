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
struct EnthalpyData
{
    double h_G = nan;
    double h_L = nan;
    double h_S = nan;

    static auto reflect()
    {
        using Self = EnthalpyData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{R::makeReflectionData("enthalpy_gas", &Self::h_G),
                          R::makeReflectionData("enthalpy_liquid", &Self::h_L),
                          R::makeReflectionData("enthalpy_solid", &Self::h_S)};
    }
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
