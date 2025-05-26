/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
