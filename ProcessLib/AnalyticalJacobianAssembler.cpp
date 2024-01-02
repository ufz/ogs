/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AnalyticalJacobianAssembler.h"

#include "LocalAssemblerInterface.h"

namespace ProcessLib
{
void AnalyticalJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, double const t, double const dt,
    std::vector<double> const& local_x, std::vector<double> const& local_x_prev,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    local_assembler.assembleWithJacobian(t, dt, local_x, local_x_prev,
                                         local_M_data, local_K_data,
                                         local_b_data, local_Jac_data);
}

void AnalyticalJacobianAssembler::assembleWithJacobianForStaggeredScheme(
    LocalAssemblerInterface& local_assembler, double const t, double const dt,
    Eigen::VectorXd const& local_x, Eigen::VectorXd const& local_x_prev,
    int const process_id, std::vector<double>& local_M_data,
    std::vector<double>& local_K_data, std::vector<double>& local_b_data,
    std::vector<double>& local_Jac_data)
{
    local_assembler.assembleWithJacobianForStaggeredScheme(
        t, dt, local_x, local_x_prev, process_id, local_M_data, local_K_data,
        local_b_data, local_Jac_data);
}

std::unique_ptr<AbstractJacobianAssembler> AnalyticalJacobianAssembler::copy()
    const
{
    return std::make_unique<AnalyticalJacobianAssembler>(*this);
}

}  // namespace ProcessLib
