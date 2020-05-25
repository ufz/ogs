/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NonlinearSolver.h"

#include "BaseLib/Logging.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "ConvergenceCriterion.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/Exceptions.h"

namespace NumLib
{
void NonlinearSolver<NonlinearSolverTag::Picard>::
    calculateNonEquilibriumInitialResiduum(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id)
{
    if (!compensate_non_equilibrium_initial_residuum_)
    {
        return;
    }

    auto& A = NumLib::GlobalMatrixProvider::provider.getMatrix(A_id_);
    auto& rhs = NumLib::GlobalVectorProvider::provider.getVector(rhs_id_);
    equation_system_->assemble(x, x_prev, process_id);
    equation_system_->getA(A);
    equation_system_->getRhs(*x_prev[process_id], rhs);

    // r_neq = A * x - rhs
    r_neq_ = &NumLib::GlobalVectorProvider::provider.getVector();
    MathLib::LinAlg::matMult(A, *x[process_id], *r_neq_);
    MathLib::LinAlg::axpy(*r_neq_, -1.0, rhs);  // res -= rhs
}

NonlinearSolverStatus NonlinearSolver<NonlinearSolverTag::Picard>::solve(
    std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev,
    std::function<void(int, std::vector<GlobalVector*> const&)> const&
        postIterationCallback,
    int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;
    auto& sys = *equation_system_;

    auto& A =
        NumLib::GlobalMatrixProvider::provider.getMatrix(A_id_);
    auto& rhs = NumLib::GlobalVectorProvider::provider.getVector(
        rhs_id_);

    std::vector<GlobalVector*> x_new{x};
    x_new[process_id] =
        &NumLib::GlobalVectorProvider::provider.getVector(x_new_id_);
    LinAlg::copy(*x[process_id], *x_new[process_id]);  // set initial guess

    bool error_norms_met = false;

    convergence_criterion_->preFirstIteration();

    int iteration = 1;
    for (; iteration <= maxiter_;
         ++iteration, convergence_criterion_->reset())
    {
        BaseLib::RunTime timer_dirichlet;
        double time_dirichlet = 0.0;

        BaseLib::RunTime time_iteration;
        time_iteration.start();

        timer_dirichlet.start();
        sys.computeKnownSolutions(*x_new[process_id], process_id);
        sys.applyKnownSolutions(*x_new[process_id]);
        time_dirichlet += timer_dirichlet.elapsed();

        sys.preIteration(iteration, *x_new[process_id]);

        BaseLib::RunTime time_assembly;
        time_assembly.start();
        sys.assemble(x_new, x_prev, process_id);
        sys.getA(A);
        sys.getRhs(*x_prev[process_id], rhs);
        INFO("[time] Assembly took {:g} s.", time_assembly.elapsed());

        // Subract non-equilibrium initial residuum if set
        if (r_neq_ != nullptr)
        {
            LinAlg::axpy(rhs, -1, *r_neq_);
        }

        timer_dirichlet.start();
        sys.applyKnownSolutionsPicard(A, rhs, *x_new[process_id]);
        time_dirichlet += timer_dirichlet.elapsed();
        INFO("[time] Applying Dirichlet BCs took {:g} s.", time_dirichlet);

        if (!sys.isLinear() && convergence_criterion_->hasResidualCheck()) {
            GlobalVector res;
            LinAlg::matMult(A, *x_new[process_id], res);  // res = A * x_new
            LinAlg::axpy(res, -1.0, rhs);   // res -= rhs
            convergence_criterion_->checkResidual(res);
        }

        BaseLib::RunTime time_linear_solver;
        time_linear_solver.start();
        bool iteration_succeeded =
            linear_solver_.solve(A, rhs, *x_new[process_id]);
        INFO("[time] Linear solver took {:g} s.", time_linear_solver.elapsed());

        if (!iteration_succeeded)
        {
            ERR("Picard: The linear solver failed.");
        }
        else
        {
            if (postIterationCallback)
            {
                postIterationCallback(iteration, x_new);
            }

            switch (sys.postIteration(*x_new[process_id]))
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
                    LinAlg::copy(*x_new[process_id], *x[process_id]);
                    break;
                case IterationResult::REPEAT_ITERATION:
                    INFO(
                        "Picard: The postIteration() hook decided that this "
                        "iteration has to be repeated.");
                    LinAlg::copy(
                        *x[process_id],
                        *x_new[process_id]);  // throw the iteration result away
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
            if (convergence_criterion_->hasDeltaXCheck()) {
                GlobalVector minus_delta_x(*x[process_id]);
                LinAlg::axpy(minus_delta_x, -1.0,
                             *x_new[process_id]);  // minus_delta_x = x - x_new
                convergence_criterion_->checkDeltaX(minus_delta_x,
                                                    *x_new[process_id]);
            }

            error_norms_met = convergence_criterion_->isSatisfied();
        }

        // Update x s.t. in the next iteration we will compute the right delta x
        LinAlg::copy(*x_new[process_id], *x[process_id]);

        INFO("[time] Iteration #{:d} took {:g} s.", iteration,
             time_iteration.elapsed());

        if (error_norms_met)
        {
            break;
        }

        // Avoid increment of the 'iteration' if the error norms are not met,
        // but maximum number of iterations is reached.
        if (iteration >= maxiter_)
        {
            break;
        }
    }

    if (iteration > maxiter_)
    {
        ERR("Picard: Could not solve the given nonlinear system within {:d} "
            "iterations",
            maxiter_);
    }

    NumLib::GlobalMatrixProvider::provider.releaseMatrix(A);
    NumLib::GlobalVectorProvider::provider.releaseVector(rhs);
    NumLib::GlobalVectorProvider::provider.releaseVector(*x_new[process_id]);

    return {error_norms_met, iteration};
}

void NonlinearSolver<NonlinearSolverTag::Newton>::
    calculateNonEquilibriumInitialResiduum(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id)
{
    if (!compensate_non_equilibrium_initial_residuum_)
    {
        return;
    }

    equation_system_->assemble(x, x_prev, process_id);
    r_neq_ = &NumLib::GlobalVectorProvider::provider.getVector();
    equation_system_->getResidual(*x[process_id], *x_prev[process_id], *r_neq_);
}

NonlinearSolverStatus NonlinearSolver<NonlinearSolverTag::Newton>::solve(
    std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev,
    std::function<void(int, std::vector<GlobalVector*> const&)> const&
        postIterationCallback,
    int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;
    auto& sys = *equation_system_;

    auto& res = NumLib::GlobalVectorProvider::provider.getVector(
        res_id_);
    auto& minus_delta_x =
        NumLib::GlobalVectorProvider::provider.getVector(
            minus_delta_x_id_);
    auto& J =
        NumLib::GlobalMatrixProvider::provider.getMatrix(J_id_);

    bool error_norms_met = false;

    // TODO be more efficient
    // init minus_delta_x to the right size
    LinAlg::copy(*x[process_id], minus_delta_x);

    convergence_criterion_->preFirstIteration();

    int iteration = 1;
    for (; iteration <= maxiter_;
         ++iteration, convergence_criterion_->reset())
    {
        BaseLib::RunTime timer_dirichlet;
        double time_dirichlet = 0.0;

        BaseLib::RunTime time_iteration;
        time_iteration.start();

        timer_dirichlet.start();
        sys.computeKnownSolutions(*x[process_id], process_id);
        sys.applyKnownSolutions(*x[process_id]);
        time_dirichlet += timer_dirichlet.elapsed();

        sys.preIteration(iteration, *x[process_id]);

        BaseLib::RunTime time_assembly;
        time_assembly.start();
        try
        {
            sys.assemble(x, x_prev, process_id);
        }
        catch (AssemblyException const& e)
        {
            ERR("Abort nonlinear iteration. Repeating timestep. Reason: {:s}",
                e.what());
            error_norms_met = false;
            iteration = maxiter_;
            break;
        }
        sys.getResidual(*x[process_id], *x_prev[process_id], res);
        sys.getJacobian(J);
        INFO("[time] Assembly took {:g} s.", time_assembly.elapsed());

        // Subract non-equilibrium initial residuum if set
        if (r_neq_ != nullptr)
            LinAlg::axpy(res, -1, *r_neq_);

        minus_delta_x.setZero();

        timer_dirichlet.start();
        sys.applyKnownSolutionsNewton(J, res, minus_delta_x);
        time_dirichlet += timer_dirichlet.elapsed();
        INFO("[time] Applying Dirichlet BCs took {:g} s.", time_dirichlet);

        if (!sys.isLinear() && convergence_criterion_->hasResidualCheck())
        {
            convergence_criterion_->checkResidual(res);
        }

        BaseLib::RunTime time_linear_solver;
        time_linear_solver.start();
        bool iteration_succeeded = linear_solver_.solve(J, res, minus_delta_x);
        INFO("[time] Linear solver took {:g} s.", time_linear_solver.elapsed());

        if (!iteration_succeeded)
        {
            ERR("Newton: The linear solver failed.");
        }
        else
        {
            // TODO could be solved in a better way
            // cf.
            // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecWAXPY.html

            // Copy pointers, replace the one for the given process id.
            std::vector<GlobalVector*> x_new{x};
            x_new[process_id] =
                &NumLib::GlobalVectorProvider::provider.getVector(
                    *x[process_id], x_new_id_);
            LinAlg::axpy(*x_new[process_id], -damping_, minus_delta_x);

            if (postIterationCallback)
            {
                postIterationCallback(iteration, x_new);
            }

            switch (sys.postIteration(*x_new[process_id]))
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
                    NumLib::GlobalVectorProvider::provider.releaseVector(
                        *x_new[process_id]);
                    continue;  // That throws the iteration result away.
            }

            LinAlg::copy(*x_new[process_id],
                         *x[process_id]);  // copy new solution to x
            NumLib::GlobalVectorProvider::provider.releaseVector(
                *x_new[process_id]);
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
            if (convergence_criterion_->hasDeltaXCheck()) {
                // Note: x contains the new solution!
                convergence_criterion_->checkDeltaX(minus_delta_x,
                                                    *x[process_id]);
            }

            error_norms_met = convergence_criterion_->isSatisfied();
        }

        INFO("[time] Iteration #{:d} took {:g} s.", iteration,
             time_iteration.elapsed());

        if (error_norms_met)
        {
            break;
        }

        // Avoid increment of the 'iteration' if the error norms are not met,
        // but maximum number of iterations is reached.
        if (iteration >= maxiter_)
        {
            break;
        }
    }

    if (iteration > maxiter_)
    {
        ERR("Newton: Could not solve the given nonlinear system within {:d} "
            "iterations",
            maxiter_);
    }

    NumLib::GlobalMatrixProvider::provider.releaseMatrix(J);
    NumLib::GlobalVectorProvider::provider.releaseVector(res);
    NumLib::GlobalVectorProvider::provider.releaseVector(
        minus_delta_x);

    return {error_norms_met, iteration};
}

std::pair<std::unique_ptr<NonlinearSolverBase>, NonlinearSolverTag>
createNonlinearSolver(GlobalLinearSolver& linear_solver,
                      BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__type}
    auto const type = config.getConfigParameter<std::string>("type");
    //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__max_iter}
    auto const max_iter = config.getConfigParameter<int>("max_iter");

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
                "{:g}.",
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
}  // namespace NumLib
