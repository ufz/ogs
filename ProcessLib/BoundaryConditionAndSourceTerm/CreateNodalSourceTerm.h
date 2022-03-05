/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
namespace ParameterLib
{
struct ParameterBase;
}
namespace ProcessLib
{
class SourceTerm;
}

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createNodalSourceTerm(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t source_term_mesh_id, const int variable_id,
    const int component_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
