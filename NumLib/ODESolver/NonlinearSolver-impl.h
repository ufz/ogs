#pragma once

#include <logog/include/logog.hpp>

// for debugging
// #include <iostream>

#include "MathLib/LinAlg/BLAS.h"

#include "NonlinearSolver.h"

#include "MathLib/LinAlg/VectorNorms.h"

// TODO: change
#include "ODETypes.h" // for one shot linear solver


namespace NumLib
{

template<typename Matrix, typename Vector>
void
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Picard>::
assemble(NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard> &sys, Vector &x)
{
    sys.assembleMatricesPicard(x);
}

template<typename Matrix, typename Vector>
bool
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Picard>::
solve(NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard> &sys, Vector &x)
{
    namespace BLAS = MathLib::BLAS;

    Matrix A; Vector rhs;
    bool success = false;

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleMatricesPicard(x);
        sys.getA(A);
        sys.getRhs(rhs);

        // std::cout << "A:\n" << Eigen::MatrixXd(A) << "\n";
        // std::cout << "rhs:\n" << rhs << "\n\n";

        oneShotLinearSolve(A, rhs, _x_new);

        BLAS::aypx(x, -1.0, _x_new); // x = _x_new - x
        auto const error = norm(x);
        // INFO("  picard iteration %u error: %e", iteration, error);

        x = _x_new;

        if (error < _tol) {
            success = true;
            break;
        }

        if (sys.isLinear()) {
            // INFO("  picard linear system. not looping");
            success = true;
            break;
        }
    }

    return success;
}


template<typename Matrix, typename Vector>
void
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Newton>::
assemble(NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton> &sys, Vector &x)
{
    sys.assembleResidualNewton(x);
}

template<typename Matrix, typename Vector>
bool
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Newton>::
solve(NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton> &sys, Vector &x)
{
    namespace BLAS = MathLib::BLAS;

    Matrix J; Vector res;
    bool success = false;

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleResidualNewton(x);
        sys.getResidual(x, res);

        // std::cout << "  res:\n" << res << std::endl;

        if (norm(res) < _tol) {
            success = true;
            break;
        }

        sys.assembleJacobian(x);
        sys.getJacobian(J);

        // std::cout << "  J:\n" << Eigen::MatrixXd(J) << std::endl;

        oneShotLinearSolve(J, res, _minus_delta_x);

        // auto const dx_norm = _minus_delta_x.norm();
        // INFO("  newton iteration %u, norm of delta x: %e", iteration, dx_norm);

        BLAS::axpy(x, -1.0, _minus_delta_x);

        if (sys.isLinear()) {
            // INFO("  newton linear system. not looping");
            success = true;
            break;
        }
    }

    return success;
}

}
