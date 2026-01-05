// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "BaseLib/TimeInterval.h"

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
struct DirichletBoundaryConditionConfig
{
    std::string parameter_name;
    BaseLib::TimeInterval time_interval;
};

DirichletBoundaryConditionConfig
parseDirichletBoundaryConditionWithinTimeIntervalConfig(
    BaseLib::ConfigTree const& config);

std::unique_ptr<BoundaryCondition>
createDirichletBoundaryConditionWithinTimeInterval(
    DirichletBoundaryConditionConfig const& config_args,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>&
        parameters);

}  // namespace ProcessLib
