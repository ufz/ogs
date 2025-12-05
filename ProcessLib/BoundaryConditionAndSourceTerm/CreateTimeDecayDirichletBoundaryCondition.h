/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-07-29 14:36:32
 */

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
