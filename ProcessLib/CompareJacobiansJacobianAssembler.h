/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
        : _asm1(std::move(asm1)),
          _asm2(std::move(asm2)),
          _abs_tol(abs_tol),
          _rel_tol(rel_tol),
          _fail_on_error(fail_on_error),
          _log_file(log_file_path)
    {
        _log_file.precision(std::numeric_limits<double>::digits10);
        _log_file << "#!/usr/bin/env python\n"
                     "import numpy as np\n"
                     "from numpy import nan\n"
                  << std::endl;
    }

    void assembleWithJacobian(LocalAssemblerInterface& local_assembler,
                              double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override;

private:
    std::unique_ptr<AbstractJacobianAssembler> _asm1;
    std::unique_ptr<AbstractJacobianAssembler> _asm2;

    // TODO change to matrix
    double const _abs_tol;
    double const _rel_tol;

    //! Whether to abort if the tolerances are exceeded.
    bool const _fail_on_error;

    //! Path where a Python script will be placed, which contains information
    //! about exceeded tolerances and assembled local matrices.
    std::ofstream _log_file;

    //! Counter used for identifying blocks in the \c _log_file. It is
    //! incremented upon each call of the assembly routine, i.e., for each
    //! element in each iteration etc.
    std::size_t _counter = 0;
};

std::unique_ptr<CompareJacobiansJacobianAssembler>
createCompareJacobiansJacobianAssembler(BaseLib::ConfigTree const& config);

}  // namespace ProcessLib
