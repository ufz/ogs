// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NonlinearSolver.h"

#include <spdlog/fmt/ranges.h>

#include <Eigen/Core>
#include <algorithm>

#include "AndersonAcceleration.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/RunTime.h"
#include "ConvergenceCriterion.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/Exceptions.h"

#ifdef USE_PETSC
#include "PETScNonlinearSolver.h"
#endif  // USE_PETSC

namespace NumLib
{
namespace
{
//! One entry of the Anderson acceleration history: the iterate \c x it was
//! taken from and the (possibly damped) step
//! \f$ f = \beta\,(g(x) - x) \f$ leading away from it.
//!
//! The two vectors are only ever appended, rotated and dropped together, so
//! keeping them in one entry rules out a desynchronization of the two buffers.
struct AndersonHistoryEntry
{
    GlobalVector* x;
    GlobalVector* f;
};

//! Returns \c entry's vectors to the provider.
void releaseAndersonHistoryEntry(AndersonHistoryEntry const& entry)
{
    NumLib::GlobalVectorProvider::provider.releaseVector(*entry.x);
    NumLib::GlobalVectorProvider::provider.releaseVector(*entry.f);
}
}  // namespace

namespace detail
{
#if !defined(USE_PETSC) && !defined(USE_LIS)
bool solvePicard(GlobalLinearSolver& linear_solver, GlobalMatrix& A,
                 GlobalVector& rhs, GlobalVector& x,
                 MathLib::LinearSolverBehaviour const linear_solver_behaviour)
{
    BaseLib::RunTime time_linear_solver;
    time_linear_solver.start();

    if (!linear_solver.compute(A, linear_solver_behaviour))
    {
        ERR("Picard: The linear solver failed in the compute() step.");
        return false;
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
                 MathLib::LinearSolverBehaviour const linear_solver_behaviour)
{
    if (linear_solver_behaviour ==
            MathLib::LinearSolverBehaviour::RECOMPUTE_AND_STORE ||
        linear_solver_behaviour == MathLib::LinearSolverBehaviour::REUSE)
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
    auto selected_global_indices =
        _equation_system->getIndicesOfResiduumWithoutInitialCompensation();

#ifdef USE_PETSC
    // Ghost entry with global index 0 is encoded as -global_size
    // After abs(), it appears as global_size and must be converted back to 0
    auto const global_size = _r_neq->size();
    for (auto& idx : selected_global_indices)
    {
        if (idx == global_size)
        {
            idx = 0;
        }
    }
#endif

    std::vector<double> zero_entries(selected_global_indices.size(), 0.0);
    _r_neq->set(selected_global_indices, zero_entries);
    _equation_system->setReleaseNodalForces(_r_neq, process_id);

    MathLib::LinAlg::finalizeAssembly(*_r_neq);

    NumLib::GlobalMatrixProvider::provider.releaseMatrix(A);
    NumLib::GlobalVectorProvider::provider.releaseVector(rhs);
}

NonlinearSolverStatus NonlinearSolver<NonlinearSolverTag::Picard>::solve(
    std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev,
    std::function<void(int, bool, std::vector<GlobalVector*> const&)> const&
        postIterationCallback,
    int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;
    auto& sys = *_equation_system;

    if ((_damping != 1.0 ||
         _anderson_depth >= AndersonAcceleration::min_mixing_depth) &&
        sys.isLinear())
    {
        OGS_FATAL(
            "Damping (under-relaxation) and Anderson acceleration are not "
            "compatible with a linear equation system: a single Picard step "
            "already yields the exact solution, so the mixed/damped iterate "
            "would be accepted as converged but wrong. Remove the 'damping' "
            "parameter and the 'anderson' subtree for linear problems.");
    }

    auto& A = NumLib::GlobalMatrixProvider::provider.getMatrix(_A_id);
    auto& rhs = NumLib::GlobalVectorProvider::provider.getVector(_rhs_id);

    std::vector<GlobalVector*> x_new{x};
    x_new[process_id] =
        &NumLib::GlobalVectorProvider::provider.getVector(_x_new_id);
    LinAlg::copy(*x[process_id], *x_new[process_id]);  // set initial guess

    bool error_norms_met = false;

    _convergence_criterion->preFirstIteration();

    // Anderson acceleration history. _anderson_depth of 0 and 1 = plain Picard.
    //
    // Circular buffer of history entries, oldest first (size <=
    // _anderson_depth). With beta = _damping = 1 the stored step reduces to the
    // plain residual g(x) - x; for beta < 1 every stored step is scaled by
    // beta, which leaves the mixing weights theta unchanged (beta cancels, see
    // below).
    std::vector<AndersonHistoryEntry> anderson_history;
    anderson_history.reserve(_anderson_depth);

    // Gram matrix G = F^T F of the stored steps, maintained incrementally
    // across iterations (only the newest step's row/column is recomputed).
    Eigen::MatrixXd gram(_anderson_depth, _anderson_depth);

    int iteration = 1;
    for (; iteration <= _maxiter; ++iteration, _convergence_criterion->reset())
    {
        BaseLib::RunTime timer_dirichlet;
        double time_dirichlet = 0.0;

        BaseLib::RunTime time_iteration;
        time_iteration.start();

        INFO("Iteration #{:d} started.", iteration);
        timer_dirichlet.start();
        auto& x_new_process = *x_new[process_id];
        LinAlg::setLocalAccessibleVector(x_new_process);
        sys.computeKnownSolutions(x_new_process, process_id);
        sys.applyKnownSolutions(x_new_process);
        time_dirichlet += timer_dirichlet.elapsed();

        sys.preIteration(iteration, x_new_process);

        BaseLib::RunTime time_assembly;
        time_assembly.start();
        bool mpi_rank_assembly_ok = true;
        try
        {
            sys.assemble(x_new, x_prev, process_id);
        }
        catch (AssemblyException const& e)
        {
            ERR("Abort nonlinear iteration. Repeating timestep. Reason: {:s}",
                e.what());
            error_norms_met = false;
            iteration = _maxiter;
            mpi_rank_assembly_ok = false;
        }
        if (BaseLib::MPI::anyOf(!mpi_rank_assembly_ok))
        {
            break;
        }
        sys.getA(A);
        sys.getRhs(*x_prev[process_id], rhs);

        // Normalize the linear equation system, if required
        if (sys.requiresNormalization() &&
            !_linear_solver.canSolveRectangular())
        {
            sys.getAandRhsNormalized(A, rhs);
            WARN(
                "The equation system is rectangular, but the current linear "
                "solver only supports square systems. "
                "The system will be normalized, which lead to a squared "
                "condition number and potential numerical issues. "
                "It is recommended to use a solver that supports rectangular "
                "equation systems for better numerical stability.");
        }

        INFO("[time] Assembly took {:g} s.", time_assembly.elapsed());

        // Subtract non-equilibrium initial residuum if set
        if (_r_neq != nullptr)
        {
            LinAlg::axpy(rhs, -1, *_r_neq);
        }

        auto const solver_needs_to_compute = sys.linearSolverNeedsToCompute();
        bool const solver_will_compute =
            _linear_solver.willCompute(solver_needs_to_compute);

        timer_dirichlet.start();
        sys.applyKnownSolutionsPicard(
            A, rhs, x_new_process,
            solver_will_compute
                ? MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE
                : MathLib::DirichletBCApplicationMode::
                      FAST_INCOMPLETE_MATRIX_UPDATE);
        time_dirichlet += timer_dirichlet.elapsed();
        INFO("[time] Applying Dirichlet BCs took {:g} s.", time_dirichlet);

        if (!sys.isLinear() && _convergence_criterion->hasResidualCheck())
        {
            if (!solver_will_compute)
            {
                // !solver_will_compute means that the Dirichlet BC application
                // is incomplete (i.e., A not properly modified) and the
                // computed residual is wrong.
                OGS_FATAL(
                    "Logic error. The solver skips the compute step for a "
                    "non-linear equation system.");
            }
            GlobalVector res;
            LinAlg::matMult(A, x_new_process, res);  // res = A * x_new
            LinAlg::axpy(res, -1.0, rhs);            // res -= rhs
            _convergence_criterion->checkResidual(res);
        }

        bool iteration_succeeded = detail::solvePicard(
            _linear_solver, A, rhs, x_new_process, solver_needs_to_compute);

        if (iteration_succeeded)
        {
            //   x_old         = x[process_id]   (iterate entering this step)
            //   x_new_process                   (raw Picard output g(x_old))
            // beta relaxation (always active when damping != 1):
            //   x_new = x_old + beta*(g(x_old) - x_old)
            //         = (1-beta)*x_old + beta*g(x_old)
            if (_damping != 1.0)
            {
                LinAlg::scale(x_new_process, _damping);
                LinAlg::axpy(x_new_process, 1.0 - _damping, *x[process_id]);
            }
            // Anderson acceleration (active when anderson_depth > 0):
            //   additionally mixes the last anderson_depth damped steps
            //   f_i = beta*(g(x_i) - x_i) (i.e. x_new_process - x_old computed
            //   after the beta relaxation above) to find the optimal theta
            //   minimising ||sum theta_i f_i|| s.t. sum theta_i = 1, then sets
            //   x_new = sum theta_i*(x_i + f_i).
            //   When anderson_depth == 0 only the beta relaxation above
            //   applies.
            if (_anderson_depth > 0)
            {
                // Whether the circular buffer is full and the oldest entry is
                // about to be evicted (needed for the incremental Gram update).
                bool const rotated =
                    static_cast<int>(anderson_history.size()) ==
                    _anderson_depth;
                if (!rotated)
                {
                    // The id out-params are unused: the provider allocates a
                    // fresh vector on every call and never re-fetches by id.
                    std::size_t x_id = 0u;
                    std::size_t f_id = 0u;
                    anderson_history.push_back(
                        {&NumLib::GlobalVectorProvider::provider.getVector(
                             x_id),
                         &NumLib::GlobalVectorProvider::provider.getVector(
                             f_id)});
                }
                else
                {
                    // Recycle the oldest entry as the newest one.
                    std::rotate(anderson_history.begin(),
                                anderson_history.begin() + 1,
                                anderson_history.end());
                }

                auto const& newest = anderson_history.back();

                // x = x_old, f = x_new_process - x_old
                LinAlg::copy(*x[process_id], *newest.x);
                LinAlg::copy(x_new_process, *newest.f);
                LinAlg::axpy(*newest.f, -1.0, *x[process_id]);

                // Actual window size, <= anderson_depth while the buffer fills.
                int const history_size =
                    static_cast<int>(anderson_history.size());

                // Incrementally maintain the (history_size x history_size) Gram
                // matrix G = F^T F whose columns are the stored damped steps
                // f_0 ... f_{history_size-1}. All steps but the newest are
                // unchanged from the previous iteration, so only the last
                // row/column is recomputed - history_size dot products instead
                // of a full history_size*(history_size+1)/2 rebuild. On a
                // rotate the oldest entry (index 0) was evicted, so the cached
                // block is first shifted up-left by one.
                if (rotated)
                {
                    gram.topLeftCorner(history_size - 1, history_size - 1) =
                        gram.block(1, 1, history_size - 1, history_size - 1)
                            .eval();
                }
                for (int i = 0; i < history_size; ++i)
                {
                    double const d =
                        LinAlg::dot(*anderson_history[i].f, *newest.f);
                    gram(i, history_size - 1) = d;
                    gram(history_size - 1, i) = d;
                }

                // A single stored step needs no mixing: the sum-to-one
                // constraint forces theta = (1), which just reproduces the
                // damped step already held in x_new_process.
                if (history_size >= 2)
                {
                    // Solve G theta = e (least-squares) with the constraint
                    // sum theta_i = 1 via a simple Lagrange formulation:
                    //
                    //   [ G  1 ] [ theta  ] = [ 0 ]
                    //   [ 1  0 ] [ lambda ]   [ 1 ]
                    //
                    // The beta factor scales G by beta^2 and cancels in theta,
                    // so the weights are identical to the undamped case. The
                    // Anderson update is then:
                    //   x_anderson = sum_i theta_i * (x_i + f_i)
                    //              = sum_i theta_i * (x_i +
                    //              beta*(g(x_i)-x_i)) = sum_i theta_i *
                    //              ((1-beta)*x_i + beta*g(x_i))
                    Eigen::MatrixXd const G =
                        gram.topLeftCorner(history_size, history_size);

                    Eigen::VectorXd const theta =
                        detail::computeAndersonWeights(G);

                    // Accumulate the Anderson mixed iterate directly into
                    // x_new_process. Its previous value is no longer needed:
                    // the newest step was already extracted from it above, and
                    // it is not aliased by any history entry (those are
                    // independent copies).
                    x_new_process.setZero();
                    for (int i = 0; i < history_size; ++i)
                    {
                        // x_new_process += theta_i * (x_i + f_i)
                        LinAlg::axpy(x_new_process, theta(i),
                                     *anderson_history[i].x);
                        LinAlg::axpy(x_new_process, theta(i),
                                     *anderson_history[i].f);
                    }

                    DBUG("Picard/Anderson: history size {:d}, theta=[{:.4g}]",
                         history_size,
                         fmt::join(theta.data(), theta.data() + history_size,
                                   ", "));
                }
            }
            // end Anderson acceleration block

            if (postIterationCallback)
            {
                postIterationCallback(iteration, error_norms_met, x_new);
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
                    // Drop the just-added (now stale) history entry, since we
                    // are repeating this iteration. In the full-buffer case
                    // this is the recycled slot; releasing it by reference is
                    // safe because the provider tracks vectors by pointer, not
                    // by the (unused) id.
                    //
                    // The rotation and the Gram shift performed above are not
                    // undone, and need not be: dropping the newest entry leaves
                    // the buffer holding the remaining entries in order, and
                    // the shifted top-left block of the Gram matrix is exactly
                    // their Gram matrix. The oldest entry stays evicted, which
                    // merely shortens the sliding window by one.
                    if (!anderson_history.empty())
                    {
                        releaseAndersonHistoryEntry(anderson_history.back());
                        anderson_history.pop_back();
                    }
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

    // Release Anderson history vectors.
    for (auto const& entry : anderson_history)
    {
        releaseAndersonHistoryEntry(entry);
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
    auto selected_global_indices =
        _equation_system->getIndicesOfResiduumWithoutInitialCompensation();

#ifdef USE_PETSC
    // Ghost entry with global index 0 is encoded as -global_size
    // After abs(), it appears as global_size and must be converted back to 0
    auto const global_size = _r_neq->size();
    for (auto& idx : selected_global_indices)
    {
        if (idx == global_size)
        {
            idx = 0;
        }
    }
#endif

    std::vector<double> zero_entries(selected_global_indices.size(), 0.0);
    _r_neq->set(selected_global_indices, zero_entries);
    _equation_system->setReleaseNodalForces(_r_neq, process_id);

    MathLib::LinAlg::finalizeAssembly(*_r_neq);
}

NonlinearSolverStatus NonlinearSolver<NonlinearSolverTag::Newton>::solve(
    std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev,
    std::function<void(int, bool, std::vector<GlobalVector*> const&)> const&
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

    NewtonStepContext step_ctx{sys, x_prev, process_id};

    int iteration = 1;
#if !defined(USE_PETSC) && !defined(USE_LIS)
    int next_iteration_inv_jacobian_recompute = 1;
#endif
    for (; iteration <= _maxiter; ++iteration, _convergence_criterion->reset())
    {
        BaseLib::RunTime timer_dirichlet;
        double time_dirichlet = 0.0;

        BaseLib::RunTime time_iteration;
        INFO("Iteration #{:d} started.", iteration);
        time_iteration.start();

        timer_dirichlet.start();
        sys.computeKnownSolutions(*x[process_id], process_id);
        time_dirichlet += timer_dirichlet.elapsed();

        sys.preIteration(iteration, *x[process_id]);

        BaseLib::RunTime time_assembly;
        time_assembly.start();
        bool mpi_rank_assembly_ok = true;
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
            mpi_rank_assembly_ok = false;
        }
        if (BaseLib::MPI::anyOf(!mpi_rank_assembly_ok))
        {
            break;
        }
        sys.getResidual(*x[process_id], *x_prev[process_id], res);
        sys.getJacobian(J);
        if (_tikhonov_lambda > 0.0 && iteration >= _tikhonov_starting_iteration)
        {
            J.addToDiagonal(_tikhonov_lambda);
        }
        INFO("[time] Assembly took {:g} s.", time_assembly.elapsed());

        // Subtract non-equilibrium initial residuum if set
        if (_r_neq != nullptr)
        {
            LinAlg::axpy(res, -1, *_r_neq);
        }

        minus_delta_x.setZero();

        timer_dirichlet.start();
        sys.applyKnownSolutionsNewton(J, res, *x[process_id], minus_delta_x);
        time_dirichlet += timer_dirichlet.elapsed();
        INFO("[time] Applying Dirichlet BCs took {:g} s.", time_dirichlet);

        if (!sys.isLinear() && _convergence_criterion->hasResidualCheck())
        {
            _convergence_criterion->checkResidual(res);
        }

        BaseLib::RunTime time_linear_solver;
        time_linear_solver.start();
#if !defined(USE_PETSC) && !defined(USE_LIS)
        auto linear_solver_behaviour = MathLib::LinearSolverBehaviour::REUSE;
        if (iteration == next_iteration_inv_jacobian_recompute)
        {
            linear_solver_behaviour =
                MathLib::LinearSolverBehaviour::RECOMPUTE_AND_STORE;
            next_iteration_inv_jacobian_recompute =
                next_iteration_inv_jacobian_recompute + _recompute_jacobian;
        }
        else if (_tikhonov_lambda > 0.0 &&
                 iteration == _tikhonov_starting_iteration)
        {
            // Force a refactorization so the newly added regularization term
            // is actually used by the linear solve instead of being
            // discarded by a reused, unregularized factorization.
            linear_solver_behaviour =
                MathLib::LinearSolverBehaviour::RECOMPUTE_AND_STORE;
        }

        bool iteration_succeeded = false;
        if (!_linear_solver.compute(J, linear_solver_behaviour))
        {
            ERR("Newton: The linear solver failed in the compute() step.");
        }
        else
        {
            iteration_succeeded = _linear_solver.solve(res, minus_delta_x);
        }
#else
        bool iteration_succeeded = _linear_solver.solve(J, res, minus_delta_x);
#endif
        INFO("[time] Linear solver took {:g} s.", time_linear_solver.elapsed());

        if (!iteration_succeeded)
        {
            ERR("Newton: The linear solver failed.");
        }
        else
        {
            // TODO could be solved in a better way
            // cf.
            // https://petsc.org/release/manualpages/Vec/VecWAXPY

            // Copy pointers, replace the one for the given process id.
            std::vector<GlobalVector*> x_new{x};
            x_new[process_id] =
                &NumLib::GlobalVectorProvider::provider.getVector(
                    *x[process_id], _x_new_id);
            auto const step_result = _step_strategy->applyStep(
                *x[process_id], minus_delta_x, res, J, *x_new[process_id],
                step_ctx, iteration);

            if (step_result.step_length != 1.0)
            {
                INFO("Step length: {:g}", step_result.step_length);
            }

            if (!step_result.success)
            {
                ERR("Newton: step strategy failed.");
                iteration_succeeded = false;
            }
            else if (!step_result.x_new_is_set)
            {
                LinAlg::axpy(*x_new[process_id], -1.0, minus_delta_x);
            }

            if (postIterationCallback)
            {
                postIterationCallback(iteration, error_norms_met, x_new);
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
