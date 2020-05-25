/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <fstream>
#include <limits>
#include <memory>
#include "AbstractJacobianAssembler.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace ProcessLib
{
//! Assembles the Jacobian matrix using two different Jacobian assemblers
//! and compares the assembled local Jacobian matrices.
//!
//! If the provided tolerances are exceeded, debugging information is logged in
//! the form of a Python script.
class CompareJacobiansJacobianAssembler final : public AbstractJacobianAssembler
{
public:
    CompareJacobiansJacobianAssembler(
        std::unique_ptr<AbstractJacobianAssembler>&& asm1,
        std::unique_ptr<AbstractJacobianAssembler>&& asm2,
        double abs_tol,
        double rel_tol,
        bool fail_on_error,
        std::string const& log_file_path)
        : asm1_(std::move(asm1)),
          asm2_(std::move(asm2)),
          abs_tol_(abs_tol),
          rel_tol_(rel_tol),
          fail_on_error_(fail_on_error),
          log_file_(log_file_path)
    {
        log_file_.precision(std::numeric_limits<double>::digits10);
        log_file_ << "#!/usr/bin/env python\n"
                     "import numpy as np\n"
                     "from numpy import nan\n"
                  << std::endl;
    }

    void assembleWithJacobian(LocalAssemblerInterface& local_assembler,
                              double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double dxdot_dx, const double dx_dx,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override;

private:
    std::unique_ptr<AbstractJacobianAssembler> asm1_;
    std::unique_ptr<AbstractJacobianAssembler> asm2_;

    // TODO change to matrix
    double const abs_tol_;
    double const rel_tol_;

    //! Whether to abort if the tolerances are exceeded.
    bool const fail_on_error_;

    //! Path where a Python script will be placed, which contains information
    //! about exceeded tolerances and assembled local matrices.
    std::ofstream log_file_;

    //! Counter used for identifying blocks in the \c log_file_. It is
    //! incremented upon each call of the assembly routine, i.e., for each
    //! element in each iteration etc.
    std::size_t counter_ = 0;
};

std::unique_ptr<CompareJacobiansJacobianAssembler>
createCompareJacobiansJacobianAssembler(BaseLib::ConfigTree const& config);

}  // namespace ProcessLib
