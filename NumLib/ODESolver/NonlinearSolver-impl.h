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
#include "MathLib/LinAlg/GlobalMatrixProviders.h"

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

    auto& A     = MathLib::GlobalMatrixProvider<Matrix>::provider.getMatrix(_A_id);
    auto& rhs   = MathLib::GlobalVectorProvider<Vector>::provider.getVector(_rhs_id);
    auto& x_new = MathLib::GlobalVectorProvider<Vector>::provider.getVector(_x_new_id);

    bool success = false;

    BLAS::copy(x, x_new); // set initial guess, TODO save the copy

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleMatricesPicard(x_new);
        sys.getA(A);
        sys.getRhs(rhs);

        // Here _x_new has to be used and it has to be equal to x!
        sys.applyKnownSolutionsPicard(A, rhs, x_new);

        // std::cout << "A:\n" << Eigen::MatrixXd(A) << "\n";
        // std::cout << "rhs:\n" << rhs << "\n\n";

        if (!_linear_solver.solve(A, rhs, x_new)) {
            ERR("The linear solver failed.");
            x = x_new;
            success = false;
            break;
        }

        // x is used as delta_x in order to compute the error.
        BLAS::aypx(x, -1.0, x_new); // x = _x_new - x
        auto const error = BLAS::norm2(x);
        // INFO("  picard iteration %u error: %e", iteration, error);

        // Update x s.t. in the next iteration we will compute the right delta x
        x = x_new;

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

    MathLib::GlobalMatrixProvider<Matrix>::provider.releaseMatrix(A);
    MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(rhs);
    MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(x_new);

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

    auto& res =
            MathLib::GlobalVectorProvider<Vector>::provider.getVector(_res_id);
    auto& minus_delta_x =
            MathLib::GlobalVectorProvider<Vector>::provider.getVector(_minus_delta_x_id);
    auto& J =
            MathLib::GlobalMatrixProvider<Matrix>::provider.getMatrix(_J_id);

    bool success = false;

    // TODO be more efficient
    // init _minus_delta_x to the right size and 0.0
    BLAS::copy(x, minus_delta_x);
    minus_delta_x.setZero();

    for (unsigned iteration=1; iteration<_maxiter; ++iteration)
    {
        sys.assembleResidualNewton(x);
        sys.getResidual(x, res);

        // std::cout << "  res:\n" << res << std::endl;

        // TODO streamline that, make consistent with Picard.
        if (BLAS::norm2(res) < _tol) {
            success = true;
            break;
        }

        sys.assembleJacobian(x);
        sys.getJacobian(J);
        sys.applyKnownSolutionsNewton(J, res, minus_delta_x);

        // std::cout << "  J:\n" << Eigen::MatrixXd(J) << std::endl;

        if (!_linear_solver.solve(J, res, minus_delta_x)) {
            ERR("The linear solver failed.");
            BLAS::axpy(x, -_alpha, minus_delta_x);
            success = false;
            break;
        }

        // auto const dx_norm = _minus_delta_x.norm();
        // INFO("  newton iteration %u, norm of delta x: %e", iteration, dx_norm);

        BLAS::axpy(x, -_alpha, minus_delta_x);

        if (sys.isLinear()) {
            // INFO("  newton linear system. not looping");
            success = true;
            break;
        }
    }

    MathLib::GlobalMatrixProvider<Matrix>::provider.releaseMatrix(J);
    MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(res);
    MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(minus_delta_x);

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
