/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NonlinearSolver.h"

#include <boost/algorithm/string.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "ConvergenceCriterion.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/Exceptions.h"
#include "PETScNonlinearSolver.h"

namespace NumLib
{
namespace detail
{
#if !defined(USE_PETSC) && !defined(USE_LIS)
bool solvePicard(GlobalLinearSolver& linear_solver, GlobalMatrix& A,
                 GlobalVector& rhs, GlobalVector& x,
                 bool const compute_necessary)
{
    BaseLib::RunTime time_linear_solver;
    time_linear_solver.start();

    if (compute_necessary)
    {
        if (!linear_solver.compute(A))
        {
            ERR("Picard: The linear solver failed in the compute() step.");
            return false;
        }
    }

    bool const iteration_succeeded = linear_solver.solve(rhs, x);

    INFO("[time] Linear solver took {:g} s.", time_linear_solver.elapsed());

    if (iteration_succeeded)
    {
        return true;
    }

    ERR("Picard: The linear solver failed in the solve() step.");
    return false;
}
#else
bool solvePicard(GlobalLinearSolver& linear_solver, GlobalMatrix& A,
                 GlobalVector& rhs, GlobalVector& x,
                 bool const compute_necessary)
{
    if (!compute_necessary)
    {
        WARN(
            "The performance optimization to skip the linear solver compute() "
            "step is not implemented for PETSc or LIS linear solvers.");
    }

    BaseLib::RunTime time_linear_solver;
    time_linear_solver.start();

    bool const iteration_succeeded = linear_solver.solve(A, rhs, x);

    INFO("[time] Linear solver took {:g} s.", time_linear_solver.elapsed());

    if (iteration_succeeded)
    {
        return true;
    }

    ERR("Picard: The linear solver failed in the solve() step.");
    return false;
}
#endif
}  // namespace detail

void NonlinearSolver<NonlinearSolverTag::Picard>::
    calculateNonEquilibriumInitialResiduum(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id)
{
    if (!_compensate_non_equilibrium_initial_residuum)
    {
        return;
    }

    INFO("Calculate non-equilibrium initial residuum.");

    auto& A = NumLib::GlobalMatrixProvider::provider.getMatrix(_A_id);
    auto& rhs = NumLib::GlobalVectorProvider::provider.getVector(_rhs_id);
    _equation_system->assemble(x, x_prev, process_id);
    _equation_system->getA(A);
    _equation_system->getRhs(*x_prev[process_id], rhs);

    // r_neq = A * x - rhs
    _r_neq = &NumLib::GlobalVectorProvider::provider.getVector(_r_neq_id);
    MathLib::LinAlg::matMult(A, *x[process_id], *_r_neq);
    MathLib::LinAlg::axpy(*_r_neq, -1.0, rhs);  // res -= rhs

    // Set the values of the selected entries of _r_neq, which are associated
    // with the equations that do not need initial residual compensation, to
    // zero.
    const auto selected_global_indices =
        _equation_system->getIndicesOfResiduumWithoutInitialCompensation();

    std::vector<double> zero_entries(selected_global_indices.size(), 0.0);
    _r_neq->set(selected_global_indices, zero_entries);

    MathLib::LinAlg::finalizeAssembly(*_r_neq);

    NumLib::GlobalMatrixProvider::provider.releaseMatrix(A);
    NumLib::GlobalVectorProvider::provider.releaseVector(rhs);
}

NonlinearSolverStatus NonlinearSolver<NonlinearSolverTag::Picard>::solve(
    std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev,
    std::function<void(int, std::vector<GlobalVector*> const&)> const&
        postIterationCallback,
    int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;
    auto& sys = *_equation_system;

    auto& A = NumLib::GlobalMatrixProvider::provider.getMatrix(_A_id);
    auto& rhs = NumLib::GlobalVectorProvider::provider.getVector(_rhs_id);

    std::vector<GlobalVector*> x_new{x};
    x_new[process_id] =
        &NumLib::GlobalVectorProvider::provider.getVector(_x_new_id);
    LinAlg::copy(*x[process_id], *x_new[process_id]);  // set initial guess

    bool error_norms_met = false;

    _convergence_criterion->preFirstIteration();

    int iteration = 1;
    for (; iteration <= _maxiter; ++iteration, _convergence_criterion->reset())
    {
        BaseLib::RunTime timer_dirichlet;
        double time_dirichlet = 0.0;

        BaseLib::RunTime time_iteration;
        time_iteration.start();

        timer_dirichlet.start();
        auto& x_new_process = *x_new[process_id];
        LinAlg::setLocalAccessibleVector(x_new_process);
        sys.computeKnownSolutions(x_new_process, process_id);
        sys.applyKnownSolutions(x_new_process);
        time_dirichlet += timer_dirichlet.elapsed();

        sys.preIteration(iteration, x_new_process);

        BaseLib::RunTime time_assembly;
        time_assembly.start();
        sys.assemble(x_new, x_prev, process_id);
        sys.getA(A);
        sys.getRhs(*x_prev[process_id], rhs);
        INFO("[time] Assembly took {:g} s.", time_assembly.elapsed());

        // Subtract non-equilibrium initial residuum if set
        if (_r_neq != nullptr)
        {
            LinAlg::axpy(rhs, -1, *_r_neq);
        }

        timer_dirichlet.start();
        sys.applyKnownSolutionsPicard(A, rhs, x_new_process);
        time_dirichlet += timer_dirichlet.elapsed();
        INFO("[time] Applying Dirichlet BCs took {:g} s.", time_dirichlet);

        if (!sys.isLinear() && _convergence_criterion->hasResidualCheck())
        {
            GlobalVector res;
            LinAlg::matMult(A, x_new_process, res);  // res = A * x_new
            LinAlg::axpy(res, -1.0, rhs);            // res -= rhs
            _convergence_criterion->checkResidual(res);
        }

        bool iteration_succeeded =
            detail::solvePicard(_linear_solver, A, rhs, x_new_process,
                                sys.linearSolverNeedsToCompute());

        if (iteration_succeeded)
        {
            if (postIterationCallback)
            {
                postIterationCallback(iteration, x_new);
            }

            switch (sys.postIteration(x_new_process))
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
                    LinAlg::copy(x_new_process, *x[process_id]);
                    break;
                case IterationResult::REPEAT_ITERATION:
                    INFO(
                        "Picard: The postIteration() hook decided that this "
                        "iteration has to be repeated.");
                    LinAlg::copy(
                        *x[process_id],
                        x_new_process);  // throw the iteration result away
                    continue;
            }
        }

        if (!iteration_succeeded)
        {
            // Don't compute error norms, break here.
            error_norms_met = false;
            break;
        }

        if (sys.isLinear())
        {
            error_norms_met = true;
        }
        else
        {
            if (_convergence_criterion->hasDeltaXCheck())
            {
                GlobalVector minus_delta_x(*x[process_id]);
                LinAlg::axpy(minus_delta_x, -1.0,
                             x_new_process);  // minus_delta_x = x - x_new
                _convergence_criterion->checkDeltaX(minus_delta_x,
                                                    x_new_process);
            }

            error_norms_met = _convergence_criterion->isSatisfied();
        }

        // Update x s.t. in the next iteration we will compute the right delta x
        LinAlg::copy(x_new_process, *x[process_id]);

        INFO("[time] Iteration #{:d} took {:g} s.", iteration,
             time_iteration.elapsed());

        if (error_norms_met)
        {
            break;
        }

        // Avoid increment of the 'iteration' if the error norms are not met,
        // but maximum number of iterations is reached.
        if (iteration >= _maxiter)
        {
            break;
        }
    }

    if (iteration > _maxiter)
    {
        ERR("Picard: Could not solve the given nonlinear system within {:d} "
            "iterations",
            _maxiter);
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
    if (!_compensate_non_equilibrium_initial_residuum)
    {
        return;
    }

    INFO("Calculate non-equilibrium initial residuum.");

    _equation_system->assemble(x, x_prev, process_id);
    _r_neq = &NumLib::GlobalVectorProvider::provider.getVector(_r_neq_id);
    _equation_system->getResidual(*x[process_id], *x_prev[process_id], *_r_neq);

    // Set the values of the selected entries of _r_neq, which are associated
    // with the equations that do not need initial residual compensation, to
    // zero.
    const auto selected_global_indices =
        _equation_system->getIndicesOfResiduumWithoutInitialCompensation();
    std::vector<double> zero_entries(selected_global_indices.size(), 0.0);

    _r_neq->set(selected_global_indices, zero_entries);

    MathLib::LinAlg::finalizeAssembly(*_r_neq);
}

NonlinearSolverStatus NonlinearSolver<NonlinearSolverTag::Newton>::solve(
    std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev,
    std::function<void(int, std::vector<GlobalVector*> const&)> const&
        postIterationCallback,
    int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;
    auto& sys = *_equation_system;

    auto& res = NumLib::GlobalVectorProvider::provider.getVector(_res_id);
    auto& minus_delta_x =
        NumLib::GlobalVectorProvider::provider.getVector(_minus_delta_x_id);
    auto& J = NumLib::GlobalMatrixProvider::provider.getMatrix(_J_id);

    bool error_norms_met = false;

    // TODO be more efficient
    // init minus_delta_x to the right size
    LinAlg::copy(*x[process_id], minus_delta_x);

    _convergence_criterion->preFirstIteration();

    int iteration = 1;
    for (; iteration <= _maxiter; ++iteration, _convergence_criterion->reset())
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
            iteration = _maxiter;
            break;
        }
        sys.getResidual(*x[process_id], *x_prev[process_id], res);
        sys.getJacobian(J);
        INFO("[time] Assembly took {:g} s.", time_assembly.elapsed());

        // Subtract non-equilibrium initial residuum if set
        if (_r_neq != nullptr)
            LinAlg::axpy(res, -1, *_r_neq);

        minus_delta_x.setZero();

        timer_dirichlet.start();
        sys.applyKnownSolutionsNewton(J, res, minus_delta_x);
        time_dirichlet += timer_dirichlet.elapsed();
        INFO("[time] Applying Dirichlet BCs took {:g} s.", time_dirichlet);

        if (!sys.isLinear() && _convergence_criterion->hasResidualCheck())
        {
            _convergence_criterion->checkResidual(res);
        }

        BaseLib::RunTime time_linear_solver;
        time_linear_solver.start();
        bool iteration_succeeded = _linear_solver.solve(J, res, minus_delta_x);
        INFO("[time] Linear solver took {:g} s.", time_linear_solver.elapsed());

        if (!iteration_succeeded)
        {
            ERR("Newton: The linear solver failed.");
        }
        else
        {
            // TODO could be solved in a better way
            // cf.
            // https://www.petsc.org/release/docs/manualpages/Vec/VecWAXPY.html

            // Copy pointers, replace the one for the given process id.
            std::vector<GlobalVector*> x_new{x};
            x_new[process_id] =
                &NumLib::GlobalVectorProvider::provider.getVector(
                    *x[process_id], _x_new_id);
            LinAlg::axpy(*x_new[process_id], -_damping, minus_delta_x);

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
                        "iteration has to be repeated.");
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

        if (sys.isLinear())
        {
            error_norms_met = true;
        }
        else
        {
            if (_convergence_criterion->hasDeltaXCheck())
            {
                // Note: x contains the new solution!
                _convergence_criterion->checkDeltaX(minus_delta_x,
                                                    *x[process_id]);
            }

            error_norms_met = _convergence_criterion->isSatisfied();
        }

        INFO("[time] Iteration #{:d} took {:g} s.", iteration,
             time_iteration.elapsed());

        if (error_norms_met)
        {
            break;
        }

        // Avoid increment of the 'iteration' if the error norms are not met,
        // but maximum number of iterations is reached.
        if (iteration >= _maxiter)
        {
            break;
        }
    }

    if (iteration > _maxiter)
    {
        ERR("Newton: Could not solve the given nonlinear system within {:d} "
            "iterations",
            _maxiter);
    }

    NumLib::GlobalMatrixProvider::provider.releaseMatrix(J);
    NumLib::GlobalVectorProvider::provider.releaseVector(res);
    NumLib::GlobalVectorProvider::provider.releaseVector(minus_delta_x);

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

    if (type == "Picard")
    {
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
#ifdef USE_PETSC
    if (boost::iequals(type, "PETScSNES"))
    {
        auto prefix =
            //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__prefix}
            config.getConfigParameter<std::string>("prefix", "");
        auto const tag = NonlinearSolverTag::Newton;
        using ConcreteNLS = PETScNonlinearSolver;
        return std::make_pair(std::make_unique<ConcreteNLS>(
                                  linear_solver, max_iter, std::move(prefix)),
                              tag);
    }

#endif
    OGS_FATAL("Unsupported nonlinear solver type '{:s}'.", type.c_str());
}

NonlinearSolver<NonlinearSolverTag::Picard>::~NonlinearSolver()
{
    if (_r_neq != nullptr)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*_r_neq);
    }
}

NonlinearSolver<NonlinearSolverTag::Newton>::~NonlinearSolver()
{
    if (_r_neq != nullptr)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*_r_neq);
    }
}

}  // namespace NumLib
