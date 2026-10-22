// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <vector>

#include "DampingPolicy.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "Types.h"

namespace NumLib
{
template <NonlinearSolverTag NLTag>
class NonlinearSystem;

/// State of the current Newton iteration passed to the step strategy. It
/// gives access to the current solution and to the equation system, so that
/// globalization strategies like line search or trust region methods can
/// re-evaluate the merit at trial points without owning the equation system.
struct NewtonStepContext
{
    NonlinearSystem<NonlinearSolverTag::Newton>& sys;
    /// Current solutions of all processes; the one being solved for is
    /// x[process_id].
    std::vector<GlobalVector*> const& x;
    std::vector<GlobalVector*> const& x_prev;
    /// Residual at x[process_id] after applying the Dirichlet BCs, i.e., the
    /// right-hand side of \f$ J\,(-\Delta x) = \mathit{res} \f$. Note: the
    /// Eigen linear solver with scaling enabled scales it in place.
    GlobalVector const& res;
    int process_id = 0;
    /// Current outer Newton iteration number, starting at 1.
    int iteration = 1;
    /// Non-equilibrium initial residuum, or nullptr if not compensated.
    GlobalVector const* r_neq = nullptr;
    /// Merit value \f$ \|P\,(r - J V)\|_2 \f$ at the current iterate, i.e.,
    /// the same quantity evaluateMerit() computes at trial points, see there.
    /// Set by the nonlinear solver before applyStep() is called.
    double merit = 0;
    /// Largest deviation \f$ \|V\|_\infty \f$ of the current iterate from the
    /// Dirichlet values, \f$ V_k = x_k - g_k \f$.
    ///
    /// It is nonzero, e.g., in the first iteration of a time step with changed
    /// Dirichlet values. The Newton step then imposes the new Dirichlet values
    /// through the right-hand side and distributes their change over the
    /// domain (predictor). After a full step it vanishes up to round-off or the
    /// linear solver tolerance. The merit accounts for the deviation, but
    /// strategies may still want to treat such iterations differently, e.g.,
    /// take the full step.
    double dirichlet_deviation = 0;

    /// Evaluates the merit \f$ \|P\,(r - J V)\|_2 \f$ at the trial point
    /// \p x_trial for the process being solved for, see
    /// NonlinearSystem::getDirichletFreeResidual(). Here \f$ r \f$ is the
    /// residual at \p x_trial with the non-equilibrium initial residuum
    /// subtracted, \f$ J \f$ the Jacobian at \p x_trial, and
    /// \f$ V_k = x_{\mathrm{trial},k} - g_k \f$ at the Dirichlet DOFs.
    ///
    /// That is the residual linearized at the state fulfilling the Dirichlet
    /// BCs, restricted to the free DOFs. It is independent of the linear
    /// algebra backend, also for trial points violating the Dirichlet BCs.
    ///
    /// The Dirichlet DOFs and values stay the ones of the current iteration,
    /// and the Jacobian used by the Newton iteration is a copy, hence not
    /// affected. However, the equation system and the process (e.g.,
    /// integration point states) hold the assembly at \p x_trial afterwards,
    /// until the next iteration reassembles.
    ///
    /// \return std::nullopt if the assembly failed on any MPI rank.
    std::optional<double> evaluateMerit(GlobalVector& x_trial) const;
};

/// Result returned by a step strategy's applyStep().
struct StepResult
{
    /// False if the strategy could not produce a usable step at all
    /// (e.g. line search failed all backtracks, trust region radius
    /// collapsed to zero). Newton::solve() should treat this like a
    /// linear solver failure and break.
    bool success;

    /// Effective step length taken. 1.0 = full Newton step.
    /// Informational — used for logging and trust region radius updates.
    double step_length;

    /// True if the strategy has already updated x_new internally and
    /// Newton::solve() should NOT apply minus_delta_x again.
    /// For FixedDamping this is true; for trust region variants that
    /// compute their own direction this is also true.
    bool x_new_is_set;
};

/*! Abstract globalization and step-acceptance strategy for Newton-Raphson
 *  iterations.
 *
 *  Concrete implementations (e.g. FixedDampingStrategy) decide how to scale
 *  the raw Newton direction \f$ -\Delta x \f$ before it is applied to the
 *  current solution.  The NonlinearSolver<Newton> owns one strategy instance
 *  and calls applyStep() once per iteration.
 */
class NewtonStepStrategy
{
public:
    /// Called once before the outer Newton loop starts.
    /// The default implementation is a no-op; override if the strategy
    /// needs to cache the initial solution (e.g. for trust-region radius
    /// initialisation).
    virtual void initialize(GlobalVector const& /*x0*/) {}

    /// Given the Newton direction \p minus_delta_x at the current solution
    /// \f$ x \f$ = *ctx.x[ctx.process_id], compute and apply a step to
    /// \p x_new. The strategy is free to re-evaluate the merit at trial points
    /// via ctx.evaluateMerit().
    ///
    /// \param[in]     minus_delta_x full Newton step (J^{-1} r)
    /// \param[in]     J             Jacobian at x (may be used for curvature)
    /// \param[in,out] x_new         on input a copy of x, on output the
    ///                              accepted new solution
    /// \param[in]     ctx           current iteration state and equation
    ///                              system for re-evaluation
    virtual StepResult applyStep(GlobalVector const& minus_delta_x,
                                 GlobalMatrix const& J,
                                 GlobalVector& x_new,
                                 NewtonStepContext const& ctx) = 0;

    /// Inject the damping policy from the convergence criterion.
    /// Called by NonlinearSolver::setEquationSystem() before the first solve.
    /// A null pointer means the criterion imposes no step-length constraint.
    void setDampingPolicy(DampingPolicy const* policy)
    {
        _damping_policy = policy;
    }

    virtual ~NewtonStepStrategy() = default;

protected:
    DampingPolicy const* _damping_policy = nullptr;

    /*! Returns the damping after applying the damping policy, if one is set,
     * to the base damping. */
    double applyDampingPolicy(GlobalVector const& minus_delta_x,
                              GlobalVector const& x, double damping) const
    {
        if (_damping_policy == nullptr)
        {
            return damping;
        }
        return _damping_policy->apply(minus_delta_x, x, damping);
    }

    /*! Applies the step scaled by \p damping to \p x_new. */
    static void applyDampedStep(GlobalVector& x_new,
                                GlobalVector const& minus_delta_x,
                                double damping)
    {
        MathLib::LinAlg::axpy(x_new, -damping, minus_delta_x);
    }
};

}  // namespace NumLib
