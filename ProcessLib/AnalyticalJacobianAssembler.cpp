/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AnalyticalJacobianAssembler.h"

#include "CoupledSolutionsForStaggeredScheme.h"
#include "LocalAssemblerInterface.h"

namespace ProcessLib
{
void AnalyticalJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, double const t, double const dt,
    std::vector<double> const& local_x, std::vector<double> const& local_xdot,
    const double dxdot_dx, const double dx_dx,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    local_assembler.assembleWithJacobian(t, dt, local_x, local_xdot, dxdot_dx,
                                         dx_dx, local_M_data, local_K_data,
                                         local_b_data, local_Jac_data);
}

void AnalyticalJacobianAssembler::assembleWithJacobianForStaggeredScheme(
    LocalAssemblerInterface& local_assembler, double const t, double const dt,
    Eigen::VectorXd const& local_x, Eigen::VectorXd const& local_xdot,
    const double dxdot_dx, const double dx_dx, int const process_id,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    local_assembler.assembleWithJacobianForStaggeredScheme(
        t, dt, local_x, local_xdot, dxdot_dx, dx_dx, process_id, local_M_data,
        local_K_data, local_b_data, local_Jac_data);
}

}  // namespace ProcessLib
