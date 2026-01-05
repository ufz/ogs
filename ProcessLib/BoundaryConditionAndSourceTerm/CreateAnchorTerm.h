// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
template <int GlobalDim>
std::unique_ptr<SourceTerm> createAnchorTerm(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t source_term_mesh_id, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
