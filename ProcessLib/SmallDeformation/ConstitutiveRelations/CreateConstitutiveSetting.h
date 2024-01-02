/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "SolidMechanics.h"

namespace BaseLib
{
class ConfigTree;
}
namespace ParameterLib
{
struct ParameterBase;
struct CoordinateSystem;
}  // namespace ParameterLib

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct CreateConstitutiveSetting
{
    static std::map<int,
                    std::unique_ptr<SolidConstitutiveRelation<DisplacementDim>>>
    createSolidConstitutiveRelations(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        std::optional<ParameterLib::CoordinateSystem> const&
            local_coordinate_system,
        BaseLib::ConfigTree const& config);
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
