/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DOFTableUtil.h"
#include <cassert>

namespace NumLib
{
namespace
{

template <class CalculateNorm>
double norm(GlobalVector const& x, unsigned const global_component,
             LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh, CalculateNorm calculate_norm)
{
    // TODO that also includes ghost nodes.
    double res = 0.0;
    MeshLib::MeshSubsets const& mss = dof_table.getMeshSubsets(global_component);
    for (unsigned i=0; i<mss.size(); i++)
    {
        MeshLib::MeshSubset const& ms = mss.getMeshSubset(i);
        if (ms.getMeshID() != mesh.getID())
            continue;
        for (MeshLib::Node const* node : ms.getNodes())
        {
            auto const value =
                getNodalValue(x, mesh, dof_table, node->getID(), global_component);

            res = calculate_norm(res, value);
        }
    }

    // TODO for PETSc some global accumulation is necessary.
    return res;
}

} // anonymous namespace

double getNodalValue(GlobalVector const& x, MeshLib::Mesh const& mesh,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::size_t const node_id,
                     std::size_t const global_component_id)
{
    MeshLib::Location const l{mesh.getID(), MeshLib::MeshItemType::Node,
                              node_id};

    auto const index = dof_table.getGlobalIndex(l, global_component_id);
    assert (index != NumLib::MeshComponentMap::nop);

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


double norm(GlobalVector const& x, unsigned const global_component,
            MathLib::VecNormType norm_type,
            LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    switch(norm_type)
    {
    case MathLib::VecNormType::NORM1:
        return norm(
            x, global_component, dof_table, mesh,
            [](double res, double value) { return res + std::abs(value); });
    case MathLib::VecNormType::NORM2:
        return std::sqrt(
            norm(x, global_component, dof_table, mesh,
                 [](double res, double value) { return res + value * value; }));
    case MathLib::VecNormType::INFINITY_N:
        return norm(x, global_component, dof_table, mesh,
                    [](double res, double value) {
                        return std::max(res, std::abs(value));
                    });
    default:
        OGS_FATAL("An invalid norm type has been passed.");
    }
}

}  // namespace NumLib
