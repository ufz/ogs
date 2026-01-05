// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <memory>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MeshLib
{
class Mesh;
}  // namespace MeshLib

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // namespace NumLib

namespace ParameterLib
{
struct ParameterBase;
template <typename T>
struct Parameter;
}  // namespace ParameterLib

namespace ProcessLib
{
class BoundaryCondition;

struct TimeDecayDirichletBoundaryConditionConfig
{
    std::string parameter_name;
    double lower_limit;
};

TimeDecayDirichletBoundaryConditionConfig
parseTimeDecayDirichletBoundaryConditionConfig(
    BaseLib::ConfigTree const& config);

std::unique_ptr<BoundaryCondition> createTimeDecayDirichletBoundaryCondition(
    TimeDecayDirichletBoundaryConditionConfig const& config,
    int const variable_id, int const component_id, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
