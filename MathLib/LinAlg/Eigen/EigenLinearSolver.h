/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "EigenOption.h"

namespace MathLib
{
class EigenMatrix;
class EigenVector;

class EigenLinearSolverBase;

class EigenLinearSolver final
{
public:
    /**
     * Constructor
     * @param solver_name A name used as a prefix for command line options
     *                    if there are such options available.
     * @param option Eigen linear solver options.
     */
    explicit EigenLinearSolver(std::string const& solver_name,
                               EigenOption const& option);

    ~EigenLinearSolver();

    /**
     * copy linear solvers options
     */
    void setOption(const EigenOption& option) { option_ = option; }

    /**
     * get linear solver options
     */
    EigenOption& getOption() { return option_; }

    bool solve(EigenMatrix& A, EigenVector& b, EigenVector& x);

protected:
    EigenOption option_;
    std::unique_ptr<EigenLinearSolverBase> solver_;
    void setRestart();
};

}  // namespace MathLib
