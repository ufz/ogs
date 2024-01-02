/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "AbstractJacobianAssembler.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{

//! Assembles the Jacobian matrix using a provided "analytical" method from the
//! local assembler.
class AnalyticalJacobianAssembler final : public AbstractJacobianAssembler
{
public:
    //! Assembles the Jacobian, the matrices \f$M\f$ and \f$K\f$, and the vector
    //! \f$b\f$.
    //! In this implementation the call is only forwarded to the respective
    //! method of the given \c local_assembler.
    void assembleWithJacobian(LocalAssemblerInterface& local_assembler,
                              double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override;

    void assembleWithJacobianForStaggeredScheme(
        LocalAssemblerInterface& local_assembler, double const t,
        double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

    std::unique_ptr<AbstractJacobianAssembler> copy() const override;
};

}  // namespace ProcessLib
