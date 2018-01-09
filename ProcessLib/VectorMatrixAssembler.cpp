/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
    : _jacobian_assembler(std::move(jacobian_assembler))
{
}

void VectorMatrixAssembler::preAssemble(
    const std::size_t mesh_item_id, LocalAssemblerInterface& local_assembler,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const GlobalVector& x)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    local_assembler.preAssemble(t, local_x);
}

void VectorMatrixAssembler::assemble(
    const std::size_t mesh_item_id, LocalAssemblerInterface& local_assembler,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const GlobalVector& x, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    const CoupledSolutionsForStaggeredScheme* cpl_xs)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);

    _local_M_data.clear();
    _local_K_data.clear();
    _local_b_data.clear();

    if (cpl_xs == nullptr)
    {
        auto const local_x = x.get(indices);
        local_assembler.assemble(t, local_x, _local_M_data, _local_K_data,
                                 _local_b_data);
    }
    else
    {
        std::vector<std::reference_wrapper<const std::vector<GlobalIndexType>>>
            indices_of_all_coupled_processes;
        indices_of_all_coupled_processes.reserve(cpl_xs->coupled_xs.size());
        for (std::size_t i = 0; i < cpl_xs->coupled_xs.size(); i++)
        {
            indices_of_all_coupled_processes.emplace_back(std::ref(indices));
        }

        auto local_coupled_xs0 = getPreviousLocalSolutions(
            *cpl_xs, indices_of_all_coupled_processes);
        auto local_coupled_xs =
            getCurrentLocalSolutions(*cpl_xs, indices_of_all_coupled_processes);

        ProcessLib::LocalCoupledSolutions local_coupled_solutions(
            cpl_xs->dt, cpl_xs->process_id, std::move(local_coupled_xs0),
            std::move(local_coupled_xs));

        local_assembler.assembleForStaggeredScheme(t, _local_M_data,
                                                   _local_K_data, _local_b_data,
                                                   local_coupled_solutions);
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
}

void VectorMatrixAssembler::assembleWithJacobian(
    std::size_t const mesh_item_id, LocalAssemblerInterface& local_assembler,
    NumLib::LocalToGlobalIndexMap const& dof_table, const double t,
    GlobalVector const& x, GlobalVector const& xdot, const double dxdot_dx,
    const double dx_dx, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac, const CoupledSolutionsForStaggeredScheme* cpl_xs,
    NumLib::LocalToGlobalIndexMap const* base_dof_table)
{
    // If base_dof_table != nullptr, then the coupled processes contains the
    // mechanical process, which is alway placed in the end of the coupled
    // process and always user higher order element than other process in the
    // coupling.
    auto const indices =
        ((base_dof_table == nullptr) ||
         (cpl_xs->process_id ==
              static_cast<int>(cpl_xs->coupled_xs.size()) - 1 &&
          cpl_xs != nullptr))
            ? NumLib::getIndices(mesh_item_id, dof_table)
            : NumLib::getIndices(mesh_item_id, *base_dof_table);
    auto const local_xdot = xdot.get(indices);

    _local_M_data.clear();
    _local_K_data.clear();
    _local_b_data.clear();
    _local_Jac_data.clear();

    if (cpl_xs == nullptr)
    {
        auto const local_x = x.get(indices);
        _jacobian_assembler->assembleWithJacobian(
            local_assembler, t, local_x, local_xdot, dxdot_dx, dx_dx,
            _local_M_data, _local_K_data, _local_b_data, _local_Jac_data);
    }
    else
    {
        if (base_dof_table == nullptr)
        {
            local_assembleWithJacobianForStaggeredScheme(
                t, indices, indices, local_xdot, local_assembler, dxdot_dx,
                dx_dx, cpl_xs);
        }
        else
        {
            if (cpl_xs->process_id ==
                static_cast<int>(cpl_xs->coupled_xs.size()) - 1)
            {
                const auto base_indices =
                    NumLib::getIndices(mesh_item_id, *base_dof_table);
                local_assembleWithJacobianForStaggeredScheme(
                    t, base_indices, indices, local_xdot, local_assembler,
                    dxdot_dx, dx_dx, cpl_xs);
            }
            else
            {
                const auto full_indices =
                    NumLib::getIndices(mesh_item_id, dof_table);
                local_assembleWithJacobianForStaggeredScheme(
                    t, indices, full_indices, local_xdot, local_assembler,
                    dxdot_dx, dx_dx, cpl_xs);
            }
        }
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
}

void VectorMatrixAssembler::local_assembleWithJacobianForStaggeredScheme(
    const double t, std::vector<GlobalIndexType> const& base_indices,
    std::vector<GlobalIndexType> const& full_indices,
    std::vector<double> const& local_xdot,
    LocalAssemblerInterface& local_assembler, const double dxdot_dx,
    const double dx_dx, CoupledSolutionsForStaggeredScheme const* cpl_xs)
{
    std::vector<std::reference_wrapper<const std::vector<GlobalIndexType>>>
        indices_of_all_coupled_processes;
    indices_of_all_coupled_processes.reserve(cpl_xs->coupled_xs.size());
    for (std::size_t i = 0; i < cpl_xs->coupled_xs.size() - 1; i++)
    {
        indices_of_all_coupled_processes.emplace_back(std::ref(base_indices));
    }
    indices_of_all_coupled_processes.emplace_back(std::ref(full_indices));

    auto local_coupled_xs0 =
        getPreviousLocalSolutions(*cpl_xs, indices_of_all_coupled_processes);
    auto local_coupled_xs =
        getCurrentLocalSolutions(*cpl_xs, indices_of_all_coupled_processes);

    ProcessLib::LocalCoupledSolutions local_coupled_solutions(
        cpl_xs->dt, cpl_xs->process_id, std::move(local_coupled_xs0),
        std::move(local_coupled_xs));

    _jacobian_assembler->assembleWithJacobianForStaggeredScheme(
        local_assembler, t, local_xdot, dxdot_dx, dx_dx, _local_M_data,
        _local_K_data, _local_b_data, _local_Jac_data, local_coupled_solutions);
}

}  // namespace ProcessLib
