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

#include "BaseLib/ConfigTree.h"
#include "MathLib/LinAlg/BLAS.h"

#include "NonlinearSolver.h"

#include "MathLib/LinAlg/VectorNorms.h"

namespace NumLib
{

template<typename Matrix, typename Vector>
void
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Picard>::
assemble(Vector const& x) const
{
    _equation_system->assembleMatricesPicard(x);
}

template<typename Matrix, typename Vector>
bool
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Picard>::
solve(Vector &x)
{
    namespace BLAS = MathLib::BLAS;
    auto& sys = *_equation_system;

    bool success = false;

    BLAS::copy(x, _x_new); // set initial guess, TODO save the copy

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleMatricesPicard(_x_new);
        sys.getA(_A);
        sys.getRhs(_rhs);

        // Here _x_new has to be used and it has to be equal to x!
        sys.applyKnownSolutionsPicard(_A, _rhs, _x_new);

        // std::cout << "A:\n" << Eigen::MatrixXd(A) << "\n";
        // std::cout << "rhs:\n" << rhs << "\n\n";

        if (!_linear_solver.solve(_A, _rhs, _x_new)) {
            ERR("The linear solver failed.");
            x = _x_new;
            success = false;
            break;
        }

        // x is used as delta_x in order to compute the error.
        BLAS::aypx(x, -1.0, _x_new); // x = _x_new - x
        auto const error = BLAS::norm2(x);
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
assemble(Vector const& x) const
{
    _equation_system->assembleResidualNewton(x);
    // TODO if the equation system would be reset to nullptr after each
    //      assemble() or solve() call, the user would be forced to set the
    //      equation every time and could not forget it.
}

template<typename Matrix, typename Vector>
bool
NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Newton>::
solve(Vector &x)
{
    namespace BLAS = MathLib::BLAS;
    auto& sys = *_equation_system;

    bool success = false;

    // TODO be more efficient
    // init _minus_delta_x to the right size and 0.0
    BLAS::copy(x, _minus_delta_x);
    _minus_delta_x.setZero();

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleResidualNewton(x);
        sys.getResidual(x, _res);

        // std::cout << "  res:\n" << res << std::endl;

        // TODO streamline that, make consistent with Picard.
        if (BLAS::norm2(_res) < _tol) {
            success = true;
            break;
        }

        sys.assembleJacobian(x);
        sys.getJacobian(_J);
        sys.applyKnownSolutionsNewton(_J, _res, _minus_delta_x);

        // std::cout << "  J:\n" << Eigen::MatrixXd(J) << std::endl;

        if (!_linear_solver.solve(_J, _res, _minus_delta_x)) {
            ERR("The linear solver failed.");
            BLAS::axpy(x, -_alpha, _minus_delta_x);
            success = false;
            break;
        }

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


template<typename Matrix, typename Vector>
std::pair<
    std::unique_ptr<NonlinearSolverBase<Matrix, Vector> >,
    NonlinearSolverTag
>
createNonlinearSolver(MathLib::LinearSolver<Matrix, Vector>& linear_solver,
                      BaseLib::ConfigTree const& config)
{
    using AbstractNLS = NonlinearSolverBase<Matrix, Vector>;

    auto const type      = config.getConfParam<std::string>("type");
    auto const tol       = config.getConfParam<double>("tol");
    auto const max_iter  = config.getConfParam<unsigned>("max_iter");

    if (type == "Picard")
    {
        auto const tag = NonlinearSolverTag::Picard;
        using ConcreteNLS = NonlinearSolver<Matrix, Vector, tag>;
        return std::make_pair(std::unique_ptr<AbstractNLS>(
            new ConcreteNLS{linear_solver, tol, max_iter}), tag);
    }
    else if (type == "Newton")
    {
        auto const tag = NonlinearSolverTag::Newton;
        using ConcreteNLS = NonlinearSolver<Matrix, Vector, tag>;
        return std::make_pair(std::unique_ptr<AbstractNLS>(
            new ConcreteNLS{linear_solver, tol, max_iter}), tag);
    }
    std::abort();
}


}
