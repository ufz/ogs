/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
assemble(NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard> &sys,
         Vector const& x) const
{
    sys.assembleMatricesPicard(x);
}

template<typename Matrix, typename Vector>
bool
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Picard>::
solve(NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard> &sys, Vector &x)
{
    namespace BLAS = MathLib::BLAS;

    bool success = false;

    BLAS::copy(x, _x_new); // set initial guess, TODO save the copy

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleMatricesPicard(_x_new);
        sys.getA(_A);
        sys.getRhs(_rhs);

        // Here _x_new has to be used and it has to be equal to x!
        sys.applyKnownComponentsPicard(_A, _rhs, _x_new);

        // std::cout << "A:\n" << Eigen::MatrixXd(A) << "\n";
        // std::cout << "rhs:\n" << rhs << "\n\n";

        oneShotLinearSolve(_A, _rhs, _x_new);

        // x is used as delta_x in order to compute the error.
        BLAS::aypx(x, -1.0, _x_new); // x = _x_new - x
        auto const error = norm(x);
        // INFO("  picard iteration %u error: %e", iteration, error);

        // Update x s.t. in the next iteration we will compute the right delta x
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
assemble(NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton> &sys,
         Vector const& x) const
{
    sys.assembleResidualNewton(x);
}

template<typename Matrix, typename Vector>
bool
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Newton>::
solve(NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton> &sys, Vector &x)
{
    namespace BLAS = MathLib::BLAS;

    bool success = false;

    // TODO init _minus_delta_x to right size and 0.0

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleResidualNewton(x);
        sys.getResidual(x, _res);

        // std::cout << "  res:\n" << res << std::endl;

        // TODO streamline the, make consistent with Picard.
        if (norm(_res) < _tol) {
            success = true;
            break;
        }

        sys.assembleJacobian(x);
        sys.getJacobian(_J);
        // TODO _minus_delta_x might not have been initialized
        sys.applyKnownComponentsNewton(_J, _res, _minus_delta_x);

        // std::cout << "  J:\n" << Eigen::MatrixXd(J) << std::endl;

        oneShotLinearSolve(_J, _res, _minus_delta_x);

        // auto const dx_norm = _minus_delta_x.norm();
        // INFO("  newton iteration %u, norm of delta x: %e", iteration, dx_norm);

        BLAS::axpy(x, -_alpha, _minus_delta_x);

        if (sys.isLinear()) {
            // INFO("  newton linear system. not looping");
            success = true;
            break;
        }
    }

    return success;
}

}
