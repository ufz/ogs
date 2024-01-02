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

    /**
     * Performs the compute() step of the Eigen linear solver.
     *
     * I.e., computes the (LU) decomposition in case of a direct solver, or
     * computes the preconditioner of an iterative solver.
     */
    bool compute(EigenMatrix& A);

    /**
     * Solves the linear system for the given right-hand side \c b and initial
     * guess \c x.
     *
     * \pre compute() must have been called before. (Not necessarily for every
     * \c x and \c b separately, but for every new/changed matrix A.
     */
    bool solve(EigenVector& b, EigenVector& x);

    /// Computes and solves in a single call.
    bool solve(EigenMatrix& A, EigenVector& b, EigenVector& x);

protected:
    EigenOption option_;
    std::unique_ptr<EigenLinearSolverBase> solver_;
    void setRestart();
    void setL();
    void setS();
    void setSmoothing();
    void setAngle();
    void setResidualUpdate();
};

}  // namespace MathLib
