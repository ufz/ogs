/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CENTRALDIFFERENCESJACOBIANASSEMBLER_H
#define PROCESSLIB_CENTRALDIFFERENCESJACOBIANASSEMBLER_H

#include <memory>
#include "AbstractJacobianAssembler.h"

namespace BaseLib
{
class ConfigTree;
}  // BaseLib

namespace ProcessLib
{
class CentralDifferencesJacobianAssembler : public AbstractJacobianAssembler
{
public:
    explicit CentralDifferencesJacobianAssembler(
        std::vector<double>&& absolute_epsilons);

    void assembleWithJacobian(
        LocalAssemblerInterface& local_assembler, double const t,
        std::vector<double> const& local_x,
        std::vector<double> const& local_xdot, const double dxdot_dx,
        const double dx_dx, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

private:
    std::vector<double> const _absolute_epsilons;
    std::vector<double> _local_M_data;
    std::vector<double> _local_K_data;
    std::vector<double> _local_b_data;
    std::vector<double> _local_x_perturbed_data;
};

std::unique_ptr<CentralDifferencesJacobianAssembler>
createCentralDifferencesJacobianAssembler(BaseLib::ConfigTree const& config);

}  // ProcessLib

#endif  // PROCESSLIB_CENTRALDIFFERENCESJACOBIANASSEMBLER_H
