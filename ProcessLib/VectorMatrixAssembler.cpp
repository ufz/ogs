/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VectorMatrixAssembler.h"

#include <cassert>
#include <functional>  // for std::reference_wrapper.

#include "NumLib/DOF/DOFTableUtil.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "LocalAssemblerInterface.h"

#include "CoupledSolutionsForStaggeredScheme.h"
#include "Process.h"

namespace ProcessLib
{
VectorMatrixAssembler::VectorMatrixAssembler(
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler)
    : jacobian_assembler_(std::move(jacobian_assembler))
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
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables,
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    for (auto dof_table : dof_tables)
    {
        indices_of_processes.emplace_back(
            NumLib::getIndices(mesh_item_id, dof_table.get()));
    }

    auto const& indices = indices_of_processes[process_id];
    local_M_data_.clear();
    local_K_data_.clear();
    local_b_data_.clear();

    if (cpl_xs == nullptr)
    {
        auto const local_x = x[process_id]->get(indices);
        auto const local_xdot = xdot[process_id]->get(indices);
        local_assembler.assemble(t, dt, local_x, local_xdot, local_M_data_,
                                 local_K_data_, local_b_data_);
    }
    else
    {
        auto local_coupled_xs0 = getCoupledLocalSolutions(cpl_xs->coupled_xs_t0,
                                                          indices_of_processes);

        auto local_coupled_xs =
            getCoupledLocalSolutions(x, indices_of_processes);

        auto const local_x = MathLib::toVector(local_coupled_xs);

        ProcessLib::LocalCoupledSolutions local_coupled_solutions(
            std::move(local_coupled_xs0));

        local_assembler.assembleForStaggeredScheme(
            t, dt, local_x, process_id, local_M_data_, local_K_data_,
            local_b_data_, local_coupled_solutions);
    }

    auto const num_r_c = indices.size();
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    if (!local_M_data_.empty())
    {
        auto const local_M = MathLib::toMatrix(local_M_data_, num_r_c, num_r_c);
        M.add(r_c_indices, local_M);
    }
    if (!local_K_data_.empty())
    {
        auto const local_K = MathLib::toMatrix(local_K_data_, num_r_c, num_r_c);
        K.add(r_c_indices, local_K);
    }
    if (!local_b_data_.empty())
    {
        assert(local_b_data_.size() == num_r_c);
        b.add(indices, local_b_data_);
    }
}

void VectorMatrixAssembler::assembleWithJacobian(
    std::size_t const mesh_item_id, LocalAssemblerInterface& local_assembler,
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables,
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac, CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    for (auto dof_table : dof_tables)
    {
        indices_of_processes.emplace_back(
            NumLib::getIndices(mesh_item_id, dof_table.get()));
    }

    auto const& indices = indices_of_processes[process_id];
    auto const local_xdot = xdot.get(indices);

    local_M_data_.clear();
    local_K_data_.clear();
    local_b_data_.clear();
    local_Jac_data_.clear();

    if (cpl_xs == nullptr)
    {
        auto const local_x = x[process_id]->get(indices);
        jacobian_assembler_->assembleWithJacobian(
            local_assembler, t, dt, local_x, local_xdot, dxdot_dx, dx_dx,
            local_M_data_, local_K_data_, local_b_data_, local_Jac_data_);
    }
    else
    {
        auto local_coupled_xs0 = getCoupledLocalSolutions(cpl_xs->coupled_xs_t0,
                                                          indices_of_processes);

        auto local_coupled_xs =
            getCoupledLocalSolutions(x, indices_of_processes);

        auto const local_x = MathLib::toVector(local_coupled_xs);

        ProcessLib::LocalCoupledSolutions local_coupled_solutions(
            std::move(local_coupled_xs0));

        jacobian_assembler_->assembleWithJacobianForStaggeredScheme(
            local_assembler, t, dt, local_x, local_xdot, dxdot_dx, dx_dx,
            process_id, local_M_data_, local_K_data_, local_b_data_,
            local_Jac_data_, local_coupled_solutions);
    }

    auto const num_r_c = indices.size();
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    if (!local_M_data_.empty())
    {
        auto const local_M = MathLib::toMatrix(local_M_data_, num_r_c, num_r_c);
        M.add(r_c_indices, local_M);
    }
    if (!local_K_data_.empty())
    {
        auto const local_K = MathLib::toMatrix(local_K_data_, num_r_c, num_r_c);
        K.add(r_c_indices, local_K);
    }
    if (!local_b_data_.empty())
    {
        assert(local_b_data_.size() == num_r_c);
        b.add(indices, local_b_data_);
    }
    if (!local_Jac_data_.empty())
    {
        auto const local_Jac =
            MathLib::toMatrix(local_Jac_data_, num_r_c, num_r_c);
        Jac.add(r_c_indices, local_Jac);
    }
    else
    {
        OGS_FATAL(
            "No Jacobian has been assembled! This might be due to programming "
            "errors in the local assembler of the current process.");
    }
}

}  // namespace ProcessLib
