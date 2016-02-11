#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <logog/include/logog.hpp>

#include <iostream>

#include "BLAS.h"

#include "NonlinearSolver.h"

#include "MathLib/LinAlg/VectorNorms.h"


void
NonlinearSolver<NonlinearSolverTag::Picard>::
assemble(NonlinearSystem<NonlinearSolverTag::Picard> &sys, Vector &x)
{
    sys.assembleMatricesPicard(x);
}

void
NonlinearSolver<NonlinearSolverTag::Picard>::
solve(NonlinearSystem<NonlinearSolverTag::Picard> &sys, Vector &x)
{
    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleMatricesPicard(x);

        auto A = sys.getA();
        auto rhs = sys.getRhs();

        // Eigen::BiCGSTAB<Matrix> linear_solver;
        // Eigen::SparseLU<Matrix> linear_solver;
        // linear_solver.compute(A);
        LinearSolver linear_solver(A);

        linear_solver.solve(rhs, _x_new);

        // _x_new = linear_solver.solveWithGuess(rhs, x);
        // _x_new = linear_solver.solve(rhs);

        BLAS::aypx(x, -1.0, _x_new); // x = _x_new - x
        auto const error = norm(x);
        // INFO("  picard iteration %u error: %e", iteration, error);

        x = _x_new;

        if (error < _tol) {
            break;
        }

        if (sys.isLinear()) {
            // INFO("  picard linear system. not looping");
            break;
        }
    }
}


void
NonlinearSolver<NonlinearSolverTag::Newton>::
assemble(NonlinearSystem<NonlinearSolverTag::Newton> &sys, Vector &x)
{
    sys.assembleResidualNewton(x);
}

void
NonlinearSolver<NonlinearSolverTag::Newton>::
solve(NonlinearSystem<NonlinearSolverTag::Newton> &sys, Vector &x)
{
    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleResidualNewton(x);

        auto res = sys.getResidual(x);

        // std::cout << "  res:\n" << res << std::endl;

        if (norm(res) < _tol) break;

        sys.assembleJacobian(x);
        auto J = sys.getJacobian();

        // std::cout << "  J:\n" << Eigen::MatrixXd(J) << std::endl;

        // Eigen::BiCGSTAB<Matrix> linear_solver;
        // Eigen::SparseLU<Matrix> linear_solver;
        // linear_solver.compute(J);
        LinearSolver linear_solver(J);

        linear_solver.solve(res, _minus_delta_x);

        // auto const dx_norm = _minus_delta_x.norm();
        // INFO("  newton iteration %u, norm of delta x: %e", iteration, dx_norm);

        x -= _minus_delta_x;

        if (sys.isLinear()) {
            // INFO("  newton linear system. not looping");
            break;
        }
    }
}

