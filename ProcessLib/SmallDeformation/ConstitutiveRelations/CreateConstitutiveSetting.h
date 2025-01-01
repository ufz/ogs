/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
namespace MeshLib
{
template <typename T>
class PropertyVector;
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
                    std::shared_ptr<SolidConstitutiveRelation<DisplacementDim>>>
    createSolidConstitutiveRelations(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        std::optional<ParameterLib::CoordinateSystem> const&
            local_coordinate_system,
        MeshLib::PropertyVector<int> const* const material_ids,
        BaseLib::ConfigTree const& config);
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
