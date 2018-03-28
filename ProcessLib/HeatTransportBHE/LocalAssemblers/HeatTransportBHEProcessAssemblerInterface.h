/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <utility>
#include <vector>

#include "BaseLib/Error.h"

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
class HeatTransportBHELocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    HeatTransportBHELocalAssemblerInterface() : _dofIndex_to_localIndex{} {}
    HeatTransportBHELocalAssemblerInterface(
        std::size_t n_local_size, std::vector<unsigned> dofIndex_to_localIndex)
        : _dofIndex_to_localIndex(std::move(dofIndex_to_localIndex))
    {
    }

    void assembleWithJacobian(double const /*t*/,
                              std::vector<double> const& /*local_x_*/,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& /*local_b_data*/,
                              std::vector<double>& /*local_Jac_data*/) override
    {
        OGS_FATAL(
            "HeatTransportBHELocalAssemblerInterface::assembleWithJacobian() "
            "is not implemented");
    }

    virtual void assembleWithJacobian(double const /*t*/,
                                      Eigen::VectorXd const& /*local_T*/,
                                      Eigen::VectorXd& /*local_b*/,
                                      Eigen::MatrixXd& /*local_A*/)
    {
        OGS_FATAL(
            "HeatTransportBHELocalAssemblerInterface::assembleWithJacobian() "
            "is not implemented");
    }

private:
    std::vector<unsigned> const _dofIndex_to_localIndex;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
