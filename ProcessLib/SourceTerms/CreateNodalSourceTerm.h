/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
namespace ProcessLib
{
class SourceTerm;
struct ParameterBase;
}  // namespace ProcessLib

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createNodalSourceTerm(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    const NumLib::LocalToGlobalIndexMap& dof_table, std::size_t mesh_id,
    const int variable_id, const int component_id,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters);

}   // namespace ProcessLib
