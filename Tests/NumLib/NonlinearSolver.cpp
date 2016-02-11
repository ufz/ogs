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
    Matrix A; Vector rhs;

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleMatricesPicard(x);
        sys.getA(A);
        sys.getRhs(rhs);

        oneShotLinearSolve(A, rhs, _x_new);

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
    Matrix J; Vector res;

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleResidualNewton(x);
        sys.getResidual(x, res);

        // std::cout << "  res:\n" << res << std::endl;

        if (norm(res) < _tol) break;

        sys.assembleJacobian(x);
        sys.getJacobian(J);

        // std::cout << "  J:\n" << Eigen::MatrixXd(J) << std::endl;

        oneShotLinearSolve(J, res, _minus_delta_x);

        // auto const dx_norm = _minus_delta_x.norm();
        // INFO("  newton iteration %u, norm of delta x: %e", iteration, dx_norm);

        BLAS::axpy(x, -1.0, _minus_delta_x);

        if (sys.isLinear()) {
            // INFO("  newton linear system. not looping");
            break;
        }
    }
}

