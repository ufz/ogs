/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include "BaseLib/Error.h"

namespace ProcessLib
{
class LocalAssemblerInterface;
struct LocalCoupledSolutions;

//! Base class for Jacobian assemblers.
class AbstractJacobianAssembler
{
public:
    //! Assembles the Jacobian, the matrices \f$M\f$ and \f$K\f$, and the vector
    //! \f$b\f$.
    virtual void assembleWithJacobian(
        LocalAssemblerInterface& local_assembler, double const t,
        std::vector<double> const& local_x,
        std::vector<double> const& local_xdot, const double dxdot_dx,
        const double dx_dx, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) = 0;

    //! Assembles the Jacobian, the matrices \f$M\f$ and \f$K\f$, and the vector
    //! \f$b\f$ with coupling.
    virtual void assembleWithJacobianForStaggeredScheme(
        LocalAssemblerInterface& /*local_assembler*/, double const /*t*/,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/,
        LocalCoupledSolutions const& /*coupled_solutions*/)
    {
        // TODO make pure virtual.
        OGS_FATAL("not implemented.");
    }

    virtual ~AbstractJacobianAssembler() = default;
};

}  // namespace ProcessLib
