#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <logog/include/logog.hpp>

#include <iostream>

#include "NonlinearSolver.h"


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

        auto const& A = sys.getA();
        auto const& rhs = sys.getRhs();

        //Eigen::BiCGSTAB<Matrix> linear_solver;
        Eigen::SparseLU<Matrix> linear_solver;
        linear_solver.compute(A);

        //_x_new = linear_solver.solveWithGuess(rhs, x);
        _x_new = linear_solver.solve(rhs);

        auto const error = (_x_new - x).norm();
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

        auto const& res = sys.getResidual(x);

        // std::cout << "  res:\n" << res << std::endl;

        if (res.norm() < _tol) break;

        sys.assembleJacobian(x);
        auto const& J = sys.getJacobian();

        // std::cout << "  J:\n" << Eigen::MatrixXd(J) << std::endl;

        // Eigen::BiCGSTAB<Matrix> linear_solver;
        Eigen::SparseLU<Matrix> linear_solver;
        linear_solver.compute(J);

        _minus_delta_x = linear_solver.solve(res);

        // auto const dx_norm = _minus_delta_x.norm();
        // INFO("  newton iteration %u, norm of delta x: %e", iteration, dx_norm);

        x -= _minus_delta_x;

        if (sys.isLinear()) {
            // INFO("  newton linear system. not looping");
            break;
        }
    }
}

