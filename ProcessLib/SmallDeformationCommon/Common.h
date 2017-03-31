/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <vector>

#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <typename LocalAssemblerInterface>
void writeNodalForces(
    MeshLib::PropertyVector<double>& nodal_forces,
    std::vector<std::unique_ptr<LocalAssemblerInterface>> const&
        local_assemblers,
    NumLib::LocalToGlobalIndexMap const& local_to_global_index_map)
{
    DBUG("Compute nodal forces for small deformation process.");

    // Zero-out the output vector before averaging.
    std::fill(std::begin(nodal_forces), std::end(nodal_forces), 0);

    GlobalExecutor::executeDereferenced(
        [](const std::size_t mesh_item_id,
           LocalAssemblerInterface& local_assembler,
           const NumLib::LocalToGlobalIndexMap& dof_table,
           std::vector<double>& node_values) {
            auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
            std::vector<double> local_data;

            local_assembler.getNodalForces(local_data);

            assert(local_data.size() == indices.size());
            for (std::size_t i = 0; i < indices.size(); ++i)
                node_values[indices[i]] += local_data[i];
        },
        local_assemblers, local_to_global_index_map, nodal_forces);
}

}  // namespace SmallDeformation
}  // namespace ProcessLib
