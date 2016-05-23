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

    bool error_norms_met = false;

    BLAS::copy(x, x_new); // set initial guess, TODO save the copy

    static unsigned unique_iteration_count = 0:


    unsigned iteration=1;
    for (; iteration<=_maxiter; ++iteration)    {
        sys.preIteration(iteration, x);

        sys.assembleMatricesPicard(x);
        sys.getA(A);
        sys.getRhs(rhs);

        // Here _x_new has to be used and it has to be equal to x!
        sys.applyKnownSolutionsPicard(A, rhs, x_new);

        // std::cout << "A:\n" << Eigen::MatrixXd(A) << "\n";
        // std::cout << "rhs:\n" << rhs << "\n\n";

		//A.write("A_" + std::to_string(iteration));
		//rhs.write("rhs_" + std::to_string(iteration));

		//++unique_iteration_count;
			
		//A.write("global_A_" + std::to_string(unique_iteration_count));

        bool iteration_succeeded = _linear_solver.solve(A, rhs, x_new);

        if (!iteration_succeeded)
        {
            ERR("Picard: The linear solver failed.");
        }
        else
        {
            switch(sys.postIteration(x_new))
            {
            case IterationResult::SUCCESS:
                // Don't copy here. The old x might still be used further below.
                // Although currently it is not.
                break;
            case IterationResult::FAILURE:
                ERR("Picard: The postIteration() hook reported a non-recoverable error.");
                iteration_succeeded = false;
                // Copy new solution to x.
                // Thereby the failed solution can be used by the caller for debugging purposes.
                BLAS::copy(x_new, x);
                break;
            case IterationResult::REPEAT_ITERATION:
                INFO("Picard: The postIteration() hook decided that this iteration"
                     " has to be repeated.");
                continue; // That throws the iteration result away.
            }
        }

        if (!iteration_succeeded) {
            // Don't compute error norms, break here.
            error_norms_met = false;
            break;
        }

        auto const norm_x = BLAS::norm2(x);
        // x is used as delta_x in order to compute the error.
        BLAS::aypx(x, -1.0, x_new); // x = _x_new - x
        auto const error_dx = BLAS::norm2(x);
        INFO("Picard: Iteration #%u |dx|=%.4e, |x|=%.4e, |dx|/|x|=%.4e,"
             " tolerance(dx)=%.4e",
             iteration, error_dx, norm_x, error_dx/norm_x, _tol);

        // Update x s.t. in the next iteration we will compute the right delta x
        BLAS::copy(x_new, x);

        if (error_dx < _tol) {
            error_norms_met = true;
            break;
        }

        if (sys.isLinear()) {
            error_norms_met = true;
            break;
        }
    }

    if (iteration > _maxiter) {
        ERR("Picard: Could not solve the given nonlinear system within %u iterations",
            _maxiter);
    }

    MathLib::GlobalMatrixProvider<Matrix>::provider.releaseMatrix(A);
    MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(rhs);
    MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(x_new);

    return error_norms_met;
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

    bool error_norms_met = false;

    // TODO be more efficient
    // init _minus_delta_x to the right size and 0.0
    BLAS::copy(x, minus_delta_x);
    minus_delta_x.setZero();

    unsigned iteration=1;
    for (; iteration<_maxiter; ++iteration)
    {
        sys.preIteration(iteration, x);

        sys.assembleResidualNewton(x);
        sys.getResidual(x, res);

        sys.assembleJacobian(x);
        sys.getJacobian(J);
        sys.applyKnownSolutionsNewton(J, res, minus_delta_x);

        auto const error_res = BLAS::norm2(res);

        // std::cout << "  J:\n" << Eigen::MatrixXd(J) << std::endl;

        bool iteration_succeeded = _linear_solver.solve(J, res, minus_delta_x);

        if (!iteration_succeeded)
        {
            ERR("Newton: The linear solver failed.");
        }
        else
        {
            // TODO could be solved in a better way
            // cf. http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecWAXPY.html
            auto& x_new =
                    MathLib::GlobalVectorProvider<Vector>::provider.getVector(x, _x_new_id);
            BLAS::axpy(x_new, -_alpha, minus_delta_x);

            switch(sys.postIteration(x_new))
            {
            case IterationResult::SUCCESS:
                break;
            case IterationResult::FAILURE:
                ERR("Newton: The postIteration() hook reported a non-recoverable error.");
                iteration_succeeded = false;
                break;
            case IterationResult::REPEAT_ITERATION:
                INFO("Newton: The postIteration() hook decided that this iteration"
                     " has to be repeated.");
                // TODO introduce some onDestroy hook.
                MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(x_new);
                continue; // That throws the iteration result away.
            }

            // TODO could be done via swap. Note: that also requires swapping the ids.
            //      Same for the Picard scheme.
            BLAS::copy(x_new, x); // copy new solution to x
            MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(x_new);
        }

        if (!iteration_succeeded) {
            // Don't compute further error norms, but break here.
            error_norms_met = false;
            break;
        }

        auto const error_dx = BLAS::norm2(minus_delta_x);
        auto const norm_x   = BLAS::norm2(x);
        INFO("Newton: Iteration #%u |dx|=%.4e, |r|=%.4e, |x|=%.4e, |dx|/|x|=%.4e,"
             " tolerance(dx)=%.4e",
             iteration, error_dx, error_res, norm_x, error_dx/norm_x, _tol);

        if (error_dx < _tol) {
            error_norms_met = true;
            break;
        }

        if (sys.isLinear()) {
            error_norms_met = true;
            break;
        }
    }

    if (iteration > _maxiter) {
        ERR("Newton: Could not solve the given nonlinear system within %u iterations",
            _maxiter);
    }

    MathLib::GlobalMatrixProvider<Matrix>::provider.releaseMatrix(J);
    MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(res);
    MathLib::GlobalVectorProvider<Vector>::provider.releaseVector(minus_delta_x);

    return error_norms_met;
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
