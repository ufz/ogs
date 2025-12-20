// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "SourceTerm.h"

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
std::unique_ptr<SourceTermBase> createEmbeddedAnchor(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    MeshLib::Mesh const& bulk_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    std::size_t source_term_mesh_id, const int variable_id);

}  // namespace ProcessLib
