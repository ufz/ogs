/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalAssemblerInterface.h"

#include <cassert>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
void LocalAssemblerInterface::assemble(
    double const /*t*/,
    double const /*dt*/,
    std::vector<double> const& /*local_x*/,
    std::vector<double> const& /*local_x_prev*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/)
{
    OGS_FATAL(
        "The assemble() function is not implemented in the local assembler.");
}

void LocalAssemblerInterface::assembleForStaggeredScheme(
    double const /*t*/, double const /*dt*/, Eigen::VectorXd const& /*local_x*/,
    Eigen::VectorXd const& /*local_x_prev*/, int const /*process_id*/,
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
    std::vector<double> const& /*local_x_prev*/,
    std::vector<double>& /*local_M_data*/,
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
    Eigen::VectorXd const& /*local_x_prev*/, int const /*process_id*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL(
        "The assembleWithJacobianForStaggeredScheme() function is not "
        "implemented in the local assembler.");
}

void LocalAssemblerInterface::computeSecondaryVariable(
    std::size_t const mesh_item_id,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, int const process_id)
{
    auto const local_x = NumLib::getLocalX(mesh_item_id, dof_tables, x);

    // Todo: A more decent way is to directly pass x_prevs as done for x
    auto const indices =
        NumLib::getIndices(mesh_item_id, *dof_tables[process_id]);
    auto const local_x_prev_vec = x_prev.get(indices);
    auto const local_x_prev = MathLib::toVector(local_x_prev_vec);

    computeSecondaryVariableConcrete(t, dt, local_x, local_x_prev);
}

void LocalAssemblerInterface::setInitialConditions(
    std::size_t const mesh_item_id,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    std::vector<GlobalVector*> const& x, double const t, int const process_id)
{
    auto local_x = NumLib::getLocalX(mesh_item_id, dof_tables, x);

    setInitialConditionsConcrete(local_x, t, process_id);
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
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    int const process_id)
{
    auto const local_x = NumLib::getLocalX(mesh_item_id, dof_tables, x);
    auto const local_x_prev =
        NumLib::getLocalX(mesh_item_id, dof_tables, x_prev);

    postTimestepConcrete(local_x, local_x_prev, t, dt, process_id);
}

void LocalAssemblerInterface::postNonLinearSolver(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    int const process_id)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x[process_id]->get(indices);
    auto const local_x_prev = x_prev[process_id]->get(indices);

    postNonLinearSolverConcrete(local_x, local_x_prev, t, dt, process_id);
}

}  // namespace ProcessLib
