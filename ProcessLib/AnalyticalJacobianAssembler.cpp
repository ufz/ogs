/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
    LocalAssemblerInterface& local_assembler, double const t,
    std::vector<double> const& local_x, std::vector<double> const& local_xdot,
    const double dxdot_dx, const double dx_dx,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data,
    std::vector<double>& local_Jac_data)
{
    local_assembler.assembleWithJacobian(t, local_x, local_xdot, dxdot_dx,
                                         dx_dx, local_M_data, local_K_data,
                                         local_b_data, local_Jac_data);
}
}  // ProcessLib
