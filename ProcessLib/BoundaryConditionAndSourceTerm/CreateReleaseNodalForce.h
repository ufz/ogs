// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

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
struct BoundaryConditionConfig;

std::string parseReleaseNodalForce(BoundaryConditionConfig const& bc_config);

std::unique_ptr<BoundaryCondition> createReleaseNodalForce(
    unsigned const global_dim, int const variable_id,
    std::string const& decay_parameter_name,
    BoundaryConditionConfig const& bc_config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
