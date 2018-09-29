/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NonlinearSolver.h"

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "ConvergenceCriterion.h"

namespace NumLib
{
void NonlinearSolver<NonlinearSolverTag::Picard>::assemble(
    GlobalVector const& x) const
{
    _equation_system->assemble(x);
}

bool NonlinearSolver<NonlinearSolverTag::Picard>::solve(
    GlobalVector& x,
    std::function<void(unsigned, GlobalVector const&)> const& postIterationCallback)
{
    namespace LinAlg = MathLib::LinAlg;
    auto& sys = *_equation_system;

    auto& A =
        NumLib::GlobalMatrixProvider::provider.getMatrix(_A_id);
    auto& rhs = NumLib::GlobalVectorProvider::provider.getVector(
        _rhs_id);
    auto& x_new = NumLib::GlobalVectorProvider::provider.getVector(_x_new_id);

    bool error_norms_met = false;

    LinAlg::copy(x, x_new);  // set initial guess

    _convergence_criterion->preFirstIteration();

    unsigned iteration = 1;
    for (; iteration <= _maxiter;
         ++iteration, _convergence_criterion->reset())
    {
        BaseLib::RunTime timer_dirichlet;
        double time_dirichlet = 0.0;

        BaseLib::RunTime time_iteration;
        time_iteration.start();

        timer_dirichlet.start();
        sys.computeKnownSolutions(x_new);
        sys.applyKnownSolutions(x_new);
        time_dirichlet += timer_dirichlet.elapsed();

        sys.preIteration(iteration, x_new);

        BaseLib::RunTime time_assembly;
        time_assembly.start();
        sys.assemble(x_new);
        sys.getA(A);
        sys.getRhs(rhs);
        INFO("[time] Assembly took %g s.", time_assembly.elapsed());

        timer_dirichlet.start();
        sys.applyKnownSolutionsPicard(A, rhs, x_new);
        time_dirichlet += timer_dirichlet.elapsed();
        INFO("[time] Applying Dirichlet BCs took %g s.", time_dirichlet);

        if (!sys.isLinear() && _convergence_criterion->hasResidualCheck()) {
            GlobalVector res;
            LinAlg::matMult(A, x_new, res);  // res = A * x_new
            LinAlg::axpy(res, -1.0, rhs);   // res -= rhs
            _convergence_criterion->checkResidual(res);
        }

        BaseLib::RunTime time_linear_solver;
        time_linear_solver.start();
        bool iteration_succeeded = _linear_solver.solve(A, rhs, x_new);
        INFO("[time] Linear solver took %g s.", time_linear_solver.elapsed());

        if (!iteration_succeeded)
        {
            ERR("Picard: The linear solver failed.");
        }
        else
        {
            if (postIterationCallback)
                postIterationCallback(iteration, x_new);

            switch (sys.postIteration(x_new))
            {
                case IterationResult::SUCCESS:
                    // Don't copy here. The old x might still be used further
                    // below. Although currently it is not.
                    break;
                case IterationResult::FAILURE:
                    ERR("Picard: The postIteration() hook reported a "
                        "non-recoverable error.");
                    iteration_succeeded = false;
                    // Copy new solution to x.
                    // Thereby the failed solution can be used by the caller for
                    // debugging purposes.
                    LinAlg::copy(x_new, x);
                    break;
                case IterationResult::REPEAT_ITERATION:
                    INFO(
                        "Picard: The postIteration() hook decided that this "
                        "iteration has to be repeated.");
                    LinAlg::copy(x, x_new);  // throw the iteration result away
                    continue;
            }
        }

        if (!iteration_succeeded)
        {
            // Don't compute error norms, break here.
            error_norms_met = false;
            break;
        }

        if (sys.isLinear()) {
            error_norms_met = true;
        } else {
            if (_convergence_criterion->hasDeltaXCheck()) {
                GlobalVector minus_delta_x(x);
                LinAlg::axpy(minus_delta_x, -1.0,
                             x_new);  // minus_delta_x = x - x_new
                _convergence_criterion->checkDeltaX(minus_delta_x, x_new);
            }

            error_norms_met = _convergence_criterion->isSatisfied();
        }

        // Update x s.t. in the next iteration we will compute the right delta x
        LinAlg::copy(x_new, x);

        INFO("[time] Iteration #%u took %g s.", iteration,
             time_iteration.elapsed());

        if (error_norms_met)
            break;
    }

    if (iteration > _maxiter)
    {
        ERR("Picard: Could not solve the given nonlinear system within %u "
            "iterations",
            _maxiter);
    }

    NumLib::GlobalMatrixProvider::provider.releaseMatrix(A);
    NumLib::GlobalVectorProvider::provider.releaseVector(rhs);
    NumLib::GlobalVectorProvider::provider.releaseVector(x_new);

    return error_norms_met;
}

void NonlinearSolver<NonlinearSolverTag::Newton>::assemble(
    GlobalVector const& x) const
{
    _equation_system->assemble(x);
    // TODO if the equation system would be reset to nullptr after each
    //      assemble() or solve() call, the user would be forced to set the
    //      equation every time and could not forget it.
}

bool NonlinearSolver<NonlinearSolverTag::Newton>::solve(
    GlobalVector& x,
    std::function<void(unsigned, GlobalVector const&)> const& postIterationCallback)
{
    namespace LinAlg = MathLib::LinAlg;
    auto& sys = *_equation_system;

    auto& res = NumLib::GlobalVectorProvider::provider.getVector(
        _res_id);
    auto& minus_delta_x =
        NumLib::GlobalVectorProvider::provider.getVector(
            _minus_delta_x_id);
    auto& J =
        NumLib::GlobalMatrixProvider::provider.getMatrix(_J_id);

    bool error_norms_met = false;

    // TODO be more efficient
    // init minus_delta_x to the right size
    LinAlg::copy(x, minus_delta_x);

    _convergence_criterion->preFirstIteration();

    unsigned iteration = 1;
    for (; iteration <= _maxiter;
         ++iteration, _convergence_criterion->reset())
    {
        BaseLib::RunTime timer_dirichlet;
        double time_dirichlet = 0.0;

        BaseLib::RunTime time_iteration;
        time_iteration.start();

        timer_dirichlet.start();
        sys.computeKnownSolutions(x);
        sys.applyKnownSolutions(x);
        time_dirichlet += timer_dirichlet.elapsed();

        sys.preIteration(iteration, x);

        BaseLib::RunTime time_assembly;
        time_assembly.start();
        sys.assemble(x);
        sys.getResidual(x, res);
        sys.getJacobian(J);
        INFO("[time] Assembly took %g s.", time_assembly.elapsed());

        minus_delta_x.setZero();

        timer_dirichlet.start();
        sys.applyKnownSolutionsNewton(J, res, minus_delta_x);
        time_dirichlet += timer_dirichlet.elapsed();
        INFO("[time] Applying Dirichlet BCs took %g s.", time_dirichlet);

        if (!sys.isLinear() && _convergence_criterion->hasResidualCheck())
            _convergence_criterion->checkResidual(res);

        BaseLib::RunTime time_linear_solver;
        time_linear_solver.start();
        bool iteration_succeeded = _linear_solver.solve(J, res, minus_delta_x);
        INFO("[time] Linear solver took %g s.", time_linear_solver.elapsed());

        if (!iteration_succeeded)
        {
            ERR("Newton: The linear solver failed.");
        }
        else
        {
            // TODO could be solved in a better way
            // cf.
            // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecWAXPY.html
            auto& x_new =
                NumLib::GlobalVectorProvider::provider.getVector(x, _x_new_id);
            LinAlg::axpy(x_new, -_damping, minus_delta_x);

            if (postIterationCallback)
                postIterationCallback(iteration, x_new);

            switch(sys.postIteration(x_new))
            {
                case IterationResult::SUCCESS:
                    break;
                case IterationResult::FAILURE:
                    ERR("Newton: The postIteration() hook reported a "
                        "non-recoverable error.");
                    iteration_succeeded = false;
                    break;
                case IterationResult::REPEAT_ITERATION:
                    INFO(
                        "Newton: The postIteration() hook decided that this "
                        "iteration"
                        " has to be repeated.");
                    // TODO introduce some onDestroy hook.
                    NumLib::GlobalVectorProvider::provider.releaseVector(x_new);
                    continue;  // That throws the iteration result away.
            }

            LinAlg::copy(x_new, x);  // copy new solution to x
            NumLib::GlobalVectorProvider::provider.releaseVector(x_new);
        }

        if (!iteration_succeeded)
        {
            // Don't compute further error norms, but break here.
            error_norms_met = false;
            break;
        }

        if (sys.isLinear()) {
            error_norms_met = true;
        } else {
            if (_convergence_criterion->hasDeltaXCheck()) {
                // Note: x contains the new solution!
                _convergence_criterion->checkDeltaX(minus_delta_x, x);
            }

            error_norms_met = _convergence_criterion->isSatisfied();
        }

        INFO("[time] Iteration #%u took %g s.", iteration,
             time_iteration.elapsed());

        if (error_norms_met)
            break;
    }

    if (iteration > _maxiter)
    {
        ERR("Newton: Could not solve the given nonlinear system within %u "
            "iterations",
            _maxiter);
    }

    NumLib::GlobalMatrixProvider::provider.releaseMatrix(J);
    NumLib::GlobalVectorProvider::provider.releaseVector(res);
    NumLib::GlobalVectorProvider::provider.releaseVector(
        minus_delta_x);

    return error_norms_met;
}

std::pair<std::unique_ptr<NonlinearSolverBase>, NonlinearSolverTag>
createNonlinearSolver(GlobalLinearSolver& linear_solver,
                      BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__type}
    auto const type = config.getConfigParameter<std::string>("type");
    //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__max_iter}
    auto const max_iter = config.getConfigParameter<unsigned>("max_iter");

    if (type == "Picard") {
        auto const tag = NonlinearSolverTag::Picard;
        using ConcreteNLS = NonlinearSolver<tag>;
        return std::make_pair(
            std::make_unique<ConcreteNLS>(linear_solver, max_iter), tag);
    }
    if (type == "Newton")
    {
        //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__damping}
        auto const damping = config.getConfigParameter<double>("damping", 1.0);
        if (damping <= 0)
        {
            OGS_FATAL(
                "The damping factor for the Newon method must be positive, got "
                "%g.",
                damping);
        }
        auto const tag = NonlinearSolverTag::Newton;
        using ConcreteNLS = NonlinearSolver<tag>;
        return std::make_pair(
            std::make_unique<ConcreteNLS>(linear_solver, max_iter, damping),
            tag);
    }
    OGS_FATAL("Unsupported nonlinear solver type");
}
}
