/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalAssemblerInterface.h"
#include <cassert>
#include "NumLib/DOF/DOFTableUtil.h"

#include "CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
void LocalAssemblerInterface::assembleWithCoupledTerm(
    double const /*t*/, std::vector<double> const& /*local_x*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    LocalCoupledSolutions const& /*coupled_solutions*/)
{
    OGS_FATAL(
        "The assembleWithCoupledTerm() function is not implemented in the "
        "local assembler.");
}

void LocalAssemblerInterface::assembleWithJacobian(
    double const /*t*/, std::vector<double> const& /*local_x*/,
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

void LocalAssemblerInterface::assembleWithJacobianAndCoupling(
    double const /*t*/, std::vector<double> const& /*local_x*/,
    std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
    const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    std::vector<double>& /*local_Jac_data*/,
    LocalCoupledSolutions const& /*coupled_solutions*/)
{
    OGS_FATAL(
        "The assembleWithJacobianAndCoupling() function is not implemented in"
        " the local assembler.");
}

void LocalAssemblerInterface::computeSecondaryVariable(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, double const t,
    GlobalVector const& x,
    CoupledSolutionsForStaggeredScheme const* coupled_term)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    if (!coupled_term)
    {
        computeSecondaryVariableConcrete(t, local_x);
    }
    else
    {
        auto const local_coupled_xs =
            getCurrentLocalSolutionsOfCoupledProcesses(coupled_term->coupled_xs,
                                                       indices);
        if (!local_coupled_xs.empty())
            computeSecondaryVariableWithCoupledProcessConcrete(
                t, local_x, local_coupled_xs);
        else
            computeSecondaryVariableConcrete(t, local_x);
    }
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
    NumLib::LocalToGlobalIndexMap const& dof_table,
    GlobalVector const& x)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    postTimestepConcrete(local_x);
}

}  // namespace ProcessLib
