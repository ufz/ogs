/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VectorMatrixAssembler.h"

#include <cassert>

#include "CoupledSolutionsForStaggeredScheme.h"
#include "LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
VectorMatrixAssembler::VectorMatrixAssembler(
    AbstractJacobianAssembler& jacobian_assembler)
    : _jacobian_assembler(jacobian_assembler)
{
}

void VectorMatrixAssembler::preAssemble(
    const std::size_t mesh_item_id, LocalAssemblerInterface& local_assembler,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    double const dt, const GlobalVector& x)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    local_assembler.preAssemble(t, dt, local_x);
}

void VectorMatrixAssembler::assemble(
    const std::size_t mesh_item_id, LocalAssemblerInterface& local_assembler,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    transform(cbegin(dof_tables), cend(dof_tables),
              back_inserter(indices_of_processes),
              [&](auto const dof_table)
              { return NumLib::getIndices(mesh_item_id, *dof_table); });

    auto const& indices = indices_of_processes[process_id];
    _local_M_data.clear();
    _local_K_data.clear();
    _local_b_data.clear();

    std::size_t const number_of_processes = x.size();
    // Monolithic scheme
    if (number_of_processes == 1)
    {
        auto const local_x = x[process_id]->get(indices);
        auto const local_x_prev = x_prev[process_id]->get(indices);
        local_assembler.assemble(t, dt, local_x, local_x_prev, _local_M_data,
                                 _local_K_data, _local_b_data);
    }
    else  // Staggered scheme
    {
        auto local_coupled_xs =
            getCoupledLocalSolutions(x, indices_of_processes);
        auto const local_x = MathLib::toVector(local_coupled_xs);

        auto local_coupled_x_prevs =
            getCoupledLocalSolutions(x_prev, indices_of_processes);
        auto const local_x_prev = MathLib::toVector(local_coupled_x_prevs);

        local_assembler.assembleForStaggeredScheme(
            t, dt, local_x, local_x_prev, process_id, _local_M_data,
            _local_K_data, _local_b_data);
    }

    auto const num_r_c = indices.size();
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    if (!_local_M_data.empty())
    {
        auto const local_M = MathLib::toMatrix(_local_M_data, num_r_c, num_r_c);
        M.add(r_c_indices, local_M);
    }
    if (!_local_K_data.empty())
    {
        auto const local_K = MathLib::toMatrix(_local_K_data, num_r_c, num_r_c);
        K.add(r_c_indices, local_K);
    }
    if (!_local_b_data.empty())
    {
        assert(_local_b_data.size() == num_r_c);
        b.add(indices, _local_b_data);
    }

    _local_output(t, process_id, mesh_item_id, _local_M_data, _local_K_data,
                  _local_b_data);
}

void VectorMatrixAssembler::assembleWithJacobian(
    std::size_t const mesh_item_id, LocalAssemblerInterface& local_assembler,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    transform(cbegin(dof_tables), cend(dof_tables),
              back_inserter(indices_of_processes),
              [&](auto const dof_table)
              { return NumLib::getIndices(mesh_item_id, *dof_table); });

    auto const& indices = indices_of_processes[process_id];

    _local_M_data.clear();
    _local_K_data.clear();
    _local_b_data.clear();
    _local_Jac_data.clear();

    std::size_t const number_of_processes = x.size();
    // Monolithic scheme
    if (number_of_processes == 1)
    {
        auto const local_x = x[process_id]->get(indices);
        auto const local_x_prev = x_prev[process_id]->get(indices);
        _jacobian_assembler.assembleWithJacobian(
            local_assembler, t, dt, local_x, local_x_prev, _local_M_data,
            _local_K_data, _local_b_data, _local_Jac_data);
    }
    else  // Staggered scheme
    {
        auto local_coupled_xs =
            getCoupledLocalSolutions(x, indices_of_processes);
        auto const local_x = MathLib::toVector(local_coupled_xs);

        auto local_coupled_x_prevs =
            getCoupledLocalSolutions(x_prev, indices_of_processes);
        auto const local_x_prev = MathLib::toVector(local_coupled_x_prevs);

        _jacobian_assembler.assembleWithJacobianForStaggeredScheme(
            local_assembler, t, dt, local_x, local_x_prev, process_id,
            _local_M_data, _local_K_data, _local_b_data, _local_Jac_data);
    }

    auto const num_r_c = indices.size();
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    if (!_local_M_data.empty())
    {
        auto const local_M = MathLib::toMatrix(_local_M_data, num_r_c, num_r_c);
        M.add(r_c_indices, local_M);
    }
    if (!_local_K_data.empty())
    {
        auto const local_K = MathLib::toMatrix(_local_K_data, num_r_c, num_r_c);
        K.add(r_c_indices, local_K);
    }
    if (!_local_b_data.empty())
    {
        assert(_local_b_data.size() == num_r_c);
        b.add(indices, _local_b_data);
    }
    if (!_local_Jac_data.empty())
    {
        auto const local_Jac =
            MathLib::toMatrix(_local_Jac_data, num_r_c, num_r_c);
        Jac.add(r_c_indices, local_Jac);
    }
    else
    {
        OGS_FATAL(
            "No Jacobian has been assembled! This might be due to programming "
            "errors in the local assembler of the current process.");
    }

    _local_output(t, process_id, mesh_item_id, _local_M_data, _local_K_data,
                  _local_b_data, &_local_Jac_data);
}

}  // namespace ProcessLib
