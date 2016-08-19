/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_ANALYTICALJACOBIANASSEMBLER_H
#define PROCESSLIB_ANALYTICALJACOBIANASSEMBLER_H

#include "AbstractJacobianAssembler.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
class AnalyticalJacobianAssembler final : public AbstractJacobianAssembler
{
public:
    void assembleWithJacobian(
        LocalAssemblerInterface& local_assembler, double const t,
        std::vector<double> const& local_x,
        std::vector<double> const& local_xdot, const double dxdot_dx,
        const double dx_dx, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;
};

}  // ProcessLib

#endif  // PROCESSLIB_ANALYTICALJACOBIANASSEMBLER_H
