/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalAssemblerInterface.h"

#include <cassert>

#include "CoupledSolutionsForStaggeredScheme.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
void LocalAssemblerInterface::assemble(
    double const /*t*/,
    double const /*dt*/,
    std::vector<double> const& /*local_x*/,
    std::vector<double> const& /*local_xdot*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/)
{
    OGS_FATAL(
        "The assemble() function is not implemented in the local assembler.");
}

void LocalAssemblerInterface::assembleForStaggeredScheme(
    double const /*t*/, double const /*dt*/, Eigen::VectorXd const& /*local_x*/,
    Eigen::VectorXd const& /*local_xdot*/, int const /*process_id*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/)
{
    OGS_FATAL(
        "The assembleForStaggeredScheme() function is not implemented in the "
        "local assembler.");
}

void LocalAssemblerInterface::assembleWithJacobian(
    double const /*t*/, double const /*dt*/,
    std::vector<double> const& /*local_x*/,
    std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
    const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL(
        "The assembleWithJacobian() function is not implemented in the local "
        "assembler.");
}

void LocalAssemblerInterface::assembleWithJacobianForStaggeredScheme(
    double const /*t*/, double const /*dt*/, Eigen::VectorXd const& /*local_x*/,
    Eigen::VectorXd const& /*local_xdot*/, const double /*dxdot_dx*/,
    const double /*dx_dx*/, int const /*process_id*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL(
        "The assembleWithJacobianForStaggeredScheme() function is not implemented in"
        " the local assembler.");
}

void LocalAssemblerInterface::computeSecondaryVariable(
    std::size_t const mesh_item_id,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_dot, int const process_id)
{
    std::vector<double> local_x_vec;

    auto const n_processes = x.size();
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        auto const indices =
            NumLib::getIndices(mesh_item_id, *dof_tables[process_id]);
        assert(!indices.empty());
        auto const local_solution = x[process_id]->get(indices);
        local_x_vec.insert(std::end(local_x_vec), std::begin(local_solution),
                           std::end(local_solution));
    }
    auto const local_x = MathLib::toVector(local_x_vec);

    // Todo: A more decent way is to directly pass x_dots as done for x
    auto const indices =
        NumLib::getIndices(mesh_item_id, *dof_tables[process_id]);
    auto const local_x_dot_vec = x_dot.get(indices);
    auto const local_x_dot = MathLib::toVector(local_x_dot_vec);

    computeSecondaryVariableConcrete(t, dt, local_x, local_x_dot);
}

void LocalAssemblerInterface::setInitialConditions(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    double const t, bool const use_monolithic_scheme, int const process_id)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    setInitialConditionsConcrete(local_x, t, use_monolithic_scheme, process_id);
}

void LocalAssemblerInterface::initialize(
    std::size_t const /*mesh_item_id*/,
    NumLib::LocalToGlobalIndexMap const& /*dof_table*/)
{
    initializeConcrete();
}

void LocalAssemblerInterface::preTimestep(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    double const t, double const delta_t)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    preTimestepConcrete(local_x, t, delta_t);
}

void LocalAssemblerInterface::postTimestep(
    std::size_t const mesh_item_id,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    std::vector<GlobalVector*> const& x, double const t, double const dt)
{
    std::vector<double> local_x_vec;

    auto const n_processes = x.size();
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        auto const indices =
            NumLib::getIndices(mesh_item_id, *dof_tables[process_id]);
        assert(!indices.empty());
        auto const local_solution = x[process_id]->get(indices);
        local_x_vec.insert(std::end(local_x_vec), std::begin(local_solution),
                           std::end(local_solution));
    }
    auto const local_x = MathLib::toVector(local_x_vec);

    postTimestepConcrete(local_x, t, dt);
}

void LocalAssemblerInterface::postNonLinearSolver(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    GlobalVector const& xdot, double const t, double const dt,
    bool const use_monolithic_scheme, int const process_id)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);
    auto const local_xdot = xdot.get(indices);

    postNonLinearSolverConcrete(local_x, local_xdot, t, dt,
                                use_monolithic_scheme, process_id);
}

}  // namespace ProcessLib
