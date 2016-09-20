/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_ABSTRACTJACOBIANASSEMBLER_H
#define PROCESSLIB_ABSTRACTJACOBIANASSEMBLER_H

#include <vector>

namespace ProcessLib
{
class LocalAssemblerInterface;

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

    virtual ~AbstractJacobianAssembler() = default;
};

}  // ProcessLib

#endif  // PROCESSLIB_ABSTRACTJACOBIANASSEMBLER_H
