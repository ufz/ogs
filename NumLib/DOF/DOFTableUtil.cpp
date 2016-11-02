/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DOFTableUtil.h"
#include <cassert>

namespace NumLib
{
double getNodalValue(GlobalVector const& x, MeshLib::Mesh const& mesh,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::size_t const node_id,
                     std::size_t const global_component_id)
{
    MeshLib::Location const l{mesh.getID(), MeshLib::MeshItemType::Node,
                              node_id};

    auto const index = dof_table.getLocalIndex(
        l, global_component_id, x.getRangeBegin(), x.getRangeEnd());
    if (index == NumLib::MeshComponentMap::nop)
        return 0.0;

    return x.get(index);
}

std::vector<GlobalIndexType> getIndices(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table)
{
    assert(dof_table.size() > mesh_item_id);
    std::vector<GlobalIndexType> indices;

    // Local matrices and vectors will always be ordered by component
    // no matter what the order of the global matrix is.
    for (unsigned c = 0; c < dof_table.getNumberOfComponents(); ++c) {
        auto const& idcs = dof_table(mesh_item_id, c).rows;
        indices.reserve(indices.size() + idcs.size());
        indices.insert(indices.end(), idcs.begin(), idcs.end());
    }

    return indices;
}

NumLib::LocalToGlobalIndexMap::RowColumnIndices getRowColumnIndices(
    std::size_t const id,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<GlobalIndexType>& indices)
{
    assert(dof_table.size() > id);
    indices.clear();

    // Local matrices and vectors will always be ordered by component,
    // no matter what the order of the global matrix is.
    for (unsigned c = 0; c < dof_table.getNumberOfComponents(); ++c)
    {
        auto const& idcs = dof_table(id, c).rows;
        indices.reserve(indices.size() + idcs.size());
        indices.insert(indices.end(), idcs.begin(), idcs.end());
    }

    return NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);
}

double norm1(GlobalVector const& x, unsigned const global_component,
             LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    // TODO that also includes ghost nodes.
    std::size_t const n_nodes = mesh.getNumberOfNodes();

    double res = 0.0;

    for (std::size_t node_id = 0; node_id < n_nodes; ++node_id) {
        auto const value =
            getNodalValue(x, mesh, dof_table, node_id, global_component);

        res += std::abs(value);
    }

    // TODO for PETSc some global accumulation is necessary.
    return res;
}

double norm2(GlobalVector const& x, unsigned const global_component,
             LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    // TODO that also includes ghost nodes.
    std::size_t const n_nodes = mesh.getNumberOfNodes();

    double res = 0.0;

    for (std::size_t node_id = 0; node_id < n_nodes; ++node_id) {
        auto const value =
            getNodalValue(x, mesh, dof_table, node_id, global_component);

        res += value*value;
    }

    // TODO for PETSc some global accumulation is necessary.
    res = std::sqrt(res);
    return res;

}

double normMax(GlobalVector const& x, unsigned const global_component,
               LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    // TODO that also includes ghost nodes.
    std::size_t const n_nodes = mesh.getNumberOfNodes();

    double res = 0.0;

    for (std::size_t node_id = 0; node_id < n_nodes; ++node_id) {
        auto const value =
            getNodalValue(x, mesh, dof_table, node_id, global_component);
        auto const abs_value = std::abs(value);
        if (abs_value > res)
            res = abs_value;
    }

    // TODO for PETSc some global accumulation is necessary.
    return res;
}

double norm(GlobalVector const& x, unsigned const global_component,
            MathLib::VecNormType norm_type,
            LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    switch(norm_type)
    {
    case MathLib::VecNormType::NORM1:
        return norm1(x, global_component, dof_table, mesh);
    case MathLib::VecNormType::NORM2:
        return norm2(x, global_component, dof_table, mesh);
    case MathLib::VecNormType::INFINITY_N:
        return normMax(x, global_component, dof_table, mesh);
    default:
        OGS_FATAL("An invalid norm type has been passed.");
    }
}

}  // namespace NumLib
