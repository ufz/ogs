// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "Types.h"

namespace NumLib
{
class DampingPolicy;
template <NonlinearSolverTag NLTag>
class NonlinearSystem;

/// Context passed to the step strategy so it can evaluate the residual
/// at trial points without owning the equation system directly.
/// Currently unused by FixedDamping/DampingReduction; reserved for the
/// in-progress line-search strategy, which needs to re-evaluate F(x) at
/// trial points.
struct NewtonStepContext
{
    NonlinearSystem<NonlinearSolverTag::Newton>& sys;
    std::vector<GlobalVector*> const& x_prev;
    int process_id;
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

    /// Given the current solution \p x, the Newton direction \p minus_delta_x,
    /// the residual \p res, and the Jacobian \p J, compute and apply a step to
    /// \p x_new. The strategy is free to re-evaluate F(x) internally.
    ///
    /// \param[in]     x             current solution
    /// \param[in]     minus_delta_x full Newton step (J^{-1} r)
    /// \param[in]     res           residual at x
    /// \param[in]     J             Jacobian at x (may be used for curvature)
    /// \param[out]    x_new         accepted new solution
    /// \param[in]     ctx           equation system context for re-evaluation
    /// \param[in]     iteration     current outer iteration number
    virtual StepResult applyStep(GlobalVector const& x,
                                 GlobalVector const& minus_delta_x,
                                 GlobalVector const& res,
                                 GlobalMatrix const& J,
                                 GlobalVector& x_new,
                                 NewtonStepContext& ctx,
                                 int iteration) = 0;
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
};

}  // namespace NumLib
