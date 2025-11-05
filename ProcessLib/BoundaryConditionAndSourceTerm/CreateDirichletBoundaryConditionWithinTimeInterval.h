/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ParameterLib
{
struct ParameterBase;
}

namespace ProcessLib
{
class BoundaryCondition;
}

namespace ProcessLib
{
std::pair<std::string, BaseLib::TimeInterval>
parseDirichletBoundaryConditionWithinTimeIntervalConfig(
    BaseLib::ConfigTree const& config);

std::unique_ptr<BoundaryCondition>
createDirichletBoundaryConditionWithinTimeInterval(
    std::string const& parameter_name,
    BaseLib::TimeInterval const& time_interval, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>&
        parameters);

}  // namespace ProcessLib
