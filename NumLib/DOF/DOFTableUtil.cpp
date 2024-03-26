/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DOFTableUtil.h"

#include <algorithm>
#include <cassert>
#include <functional>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
namespace NumLib
{
namespace
{
template <typename CalculateNorm>
double norm(GlobalVector const& x, unsigned const global_component,
            LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
            CalculateNorm calculate_norm)
{
#ifdef USE_PETSC
    x.setLocalAccessibleVector();
#endif

    double res = 0.0;
    auto const& ms = dof_table.getMeshSubset(global_component);

    assert(ms.getMeshID() == mesh.getID());
    for (MeshLib::Node const* node : ms.getNodes())
    {
        auto const value = getNonGhostNodalValue(
            x, mesh, dof_table, node->getID(), global_component);
        res = calculate_norm(res, value);
    }
    return res;
}

double norm1(GlobalVector const& x, unsigned const global_component,
             LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    double res =
        norm(x, global_component, dof_table, mesh,
             [](double res, double value) { return res + std::abs(value); });

#ifdef USE_PETSC
    double global_result = 0.0;
    MPI_Allreduce(&res, &global_result, 1, MPI_DOUBLE, MPI_SUM,
                  PETSC_COMM_WORLD);
    res = global_result;
#endif
    return res;
}

double norm2(GlobalVector const& x, unsigned const global_component,
             LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    double res =
        norm(x, global_component, dof_table, mesh,
             [](double res, double value) { return res + value * value; });

#ifdef USE_PETSC
    double global_result = 0.0;
    MPI_Allreduce(&res, &global_result, 1, MPI_DOUBLE, MPI_SUM,
                  PETSC_COMM_WORLD);
    res = global_result;
#endif
    return std::sqrt(res);
}

double normInfinity(GlobalVector const& x, unsigned const global_component,
                    LocalToGlobalIndexMap const& dof_table,
                    MeshLib::Mesh const& mesh)
{
    double res = norm(x, global_component, dof_table, mesh,
                      [](double res, double value)
                      { return std::max(res, std::abs(value)); });

#ifdef USE_PETSC
    double global_result = 0.0;
    MPI_Allreduce(&res, &global_result, 1, MPI_DOUBLE, MPI_MAX,
                  PETSC_COMM_WORLD);
    res = global_result;
#endif
    return res;
}
}  // anonymous namespace

double getNonGhostNodalValue(GlobalVector const& x, MeshLib::Mesh const& mesh,
                             NumLib::LocalToGlobalIndexMap const& dof_table,
                             std::size_t const node_id,
                             std::size_t const global_component_id)
{
    MeshLib::Location const l{mesh.getID(), MeshLib::MeshItemType::Node,
                              node_id};

    auto const index = dof_table.getGlobalIndex(l, global_component_id);
    assert(index != NumLib::MeshComponentMap::nop);

    if (index < 0)
    {  // ghost node value
        return 0.0;
    }

    return x.get(index);
}

double getNodalValue(GlobalVector const& x, MeshLib::Mesh const& mesh,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::size_t const node_id,
                     std::size_t const global_component_id)
{
    MeshLib::Location const l{mesh.getID(), MeshLib::MeshItemType::Node,
                              node_id};

    auto const index = dof_table.getGlobalIndex(l, global_component_id);
    assert(index != NumLib::MeshComponentMap::nop);

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
    for (int c = 0; c < dof_table.getNumberOfGlobalComponents(); ++c)
    {
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
    for (int c = 0; c < dof_table.getNumberOfGlobalComponents(); ++c)
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
    switch (norm_type)
    {
        case MathLib::VecNormType::NORM1:
            return norm1(x, global_component, dof_table, mesh);
        case MathLib::VecNormType::NORM2:
            return norm2(x, global_component, dof_table, mesh);
        case MathLib::VecNormType::INFINITY_N:
            return normInfinity(x, global_component, dof_table, mesh);
        default:
            OGS_FATAL("An invalid norm type has been passed.");
    }
}

std::vector<NumLib::LocalToGlobalIndexMap const*> getDOFTables(
    int const number_of_processes,
    std::function<NumLib::LocalToGlobalIndexMap const&(const int)>
        get_single_dof_table)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(number_of_processes);
    std::generate_n(std::back_inserter(dof_tables), number_of_processes,
                    [&]() { return &get_single_dof_table(dof_tables.size()); });
    return dof_tables;
}

Eigen::VectorXd getLocalX(
    std::size_t const mesh_item_id,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    std::vector<GlobalVector*> const& x)
{
    Eigen::VectorXd local_x_vec;

    auto const n_processes = x.size();
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        auto const indices =
            NumLib::getIndices(mesh_item_id, *dof_tables[process_id]);
        assert(!indices.empty());
        auto const last = local_x_vec.size();
        local_x_vec.conservativeResize(last + indices.size());
        auto const local_solution = x[process_id]->get(indices);
        assert(indices.size() == local_solution.size());
        local_x_vec.tail(local_solution.size()).noalias() =
            MathLib::toVector(local_solution);
    }
    return local_x_vec;
}

}  // namespace NumLib
