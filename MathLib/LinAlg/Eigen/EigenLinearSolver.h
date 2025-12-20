// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "EigenOption.h"
#include "MathLib/LinAlg/LinearSolverBehaviour.h"

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
    bool compute(EigenMatrix& A,
                 MathLib::LinearSolverBehaviour const linear_solver_behaviour);
    /**
     * Solves the linear system for the given right-hand side \c b and initial
     * guess \c x.
     *
     * \pre compute() must have been called before. (Not necessarily for every
     * \c x and \c b separately, but for every new/changed matrix A.
     */
    bool solve(EigenVector& b, EigenVector& x);

    /// Computes and solves in a single call.
    bool solve(EigenMatrix& A,
               EigenVector& b,
               EigenVector& x,
               MathLib::LinearSolverBehaviour const linear_solver_behaviour =
                   MathLib::LinearSolverBehaviour::RECOMPUTE);

    /// Get, if the solver can handle rectangular equation systems
    bool canSolveRectangular() const { return can_solve_rectangular_; }

    /// Tells if the solver will perform the compute step the next time it is
    /// called or if it can reuse the results from the preceding call.
    bool willCompute(
        MathLib::LinearSolverBehaviour const linear_solver_behaviour) const;

protected:
    EigenOption option_;
    std::unique_ptr<EigenLinearSolverBase> solver_;
    bool can_solve_rectangular_ = false;
    void setRestart();
    void setL();
    void setS();
    void setSmoothing();
    void setAngle();
    void setResidualUpdate();
};

}  // namespace MathLib
