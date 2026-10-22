// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "CreateTestLinearSolver.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/Exceptions.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/NumericsConfig.h"
#include "NumLib/ODESolver/ConvergenceCriterion.h"
#include "NumLib/ODESolver/NewtonStepStrategy.h"
#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "ODEs.h"
#ifdef USE_PETSC
#include "Tests/MathLib/PETScSparsityPatternUtils.h"
#endif

namespace
{
using BCIndex = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;

/*! ODE1 with Dirichlet BCs, i.e., known solutions, which can be changed
 * between calls of computeKnownSolutions().
 *
 * It is reported as nonlinear, s.t. the Newton solver iterates and checks the
 * residual.
 */
class ODEWithDirichletBCs final : public ODE1
{
public:
    std::vector<NumLib::IndexValueVector<BCIndex>> const* getKnownSolutions(
        double const /*t*/, GlobalVector const& /*x*/,
        int const /*process_id*/) const override
    {
        if (known_solutions_per_call.empty())
        {
            return &known_solutions;
        }
        auto const call = std::min(number_of_known_solutions_calls++,
                                   known_solutions_per_call.size() - 1);
        return &known_solutions_per_call[call];
    }

    void assembleWithJacobian(const double t, double const dt,
                              std::vector<GlobalVector*> const& x_curr,
                              std::vector<GlobalVector*> const& x_prev,
                              int const process_id, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        if (assembly_fails)
        {
            throw NumLib::AssemblyException("ODEWithDirichletBCs: by design.");
        }
        ODE1::assembleWithJacobian(t, dt, x_curr, x_prev, process_id, b, Jac);
    }

    bool isLinear() const override { return false; }

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int /*process_id*/) const override
    {
        return {N, N, nullptr, sparsity_pattern};
    }

    std::vector<NumLib::IndexValueVector<BCIndex>> known_solutions;
    //! If not empty, the n-th call of getKnownSolutions() returns the n-th
    //! entry instead of known_solutions, and the last one after the end.
    std::vector<std::vector<NumLib::IndexValueVector<BCIndex>>>
        known_solutions_per_call;
    mutable std::size_t number_of_known_solutions_calls = 0;
    bool assembly_fails = false;
    //! Required for PETSc matrices only.
    GlobalSparsityPattern const* sparsity_pattern = nullptr;
};

//! Records the residuals the nonlinear solver checks, and is never satisfied,
//! s.t. the solver runs for its maximum number of iterations.
class RecordingResidualCriterion final : public NumLib::ConvergenceCriterion
{
public:
    RecordingResidualCriterion()
        : NumLib::ConvergenceCriterion(MathLib::VecNormType::NORM2)
    {
    }

    bool hasDeltaXCheck() const override { return false; }
    bool hasResidualCheck() const override { return true; }

    void checkDeltaX(GlobalVector const& /*minus_delta_x*/,
                     GlobalVector const& /*x*/) override
    {
    }

    void checkResidual(GlobalVector const& residual) override
    {
        residuals.emplace_back(residual);
        _satisfied = false;
    }

    std::vector<GlobalVector> residuals;
};

//! Takes full Newton steps and records the merit values and the deviations
//! from the Dirichlet values the solver passes in each iteration.
class RecordingStepStrategy final : public NumLib::NewtonStepStrategy
{
public:
    NumLib::StepResult applyStep(GlobalVector const& /*x*/,
                                 GlobalVector const& /*minus_delta_x*/,
                                 GlobalVector const& /*res*/,
                                 GlobalMatrix const& /*J*/,
                                 GlobalVector& /*x_new*/,
                                 NumLib::NewtonStepContext& ctx,
                                 int const /*iteration*/) override
    {
        merits.push_back(ctx.merit);
        dirichlet_deviations.push_back(ctx.dirichlet_deviation);
        return {.success = true, .step_length = 1.0, .x_new_is_set = false};
    }

    std::vector<double> merits;
    std::vector<double> dirichlet_deviations;
};

using DirichletNewtonODESystem =
    NumLib::TimeDiscretizedODESystem<ODEWithDirichletBCs::ODETag,
                                     NumLib::NonlinearSolverTag::Newton>;

GlobalVector makeDirichletTestVector(std::initializer_list<double> const values)
{
    GlobalVector v(static_cast<BCIndex>(values.size()));
    MathLib::setVector(v, values);
    MathLib::LinAlg::finalizeAssembly(v);
    return v;
}
}  // namespace

// Disabled for PETSc together with the other ODE system tests, see issue
// #1989: the ODE systems used for testing do not provide a sparsity pattern,
// which PETSc requires for matrix preallocation.

/*! The initial guess violates the Dirichlet BC. In the first iteration the
 * step strategy sees the deviation, the residual at the Dirichlet DOF is
 * nonzero, and the full Newton step imposes the BC. Afterwards the deviation
 * and the residual at the Dirichlet DOF vanish up to round-off.
 */
#ifndef USE_PETSC
TEST(NumLibDirichletResidual, FullStepAfterDirichletPredictorFulfillsBC)
#else
TEST(NumLibDirichletResidual,
     DISABLED_FullStepAfterDirichletPredictorFulfillsBC)
#endif
{
    int const process_id = 0;
    double const t0 = 0.0;
    double const dt = 1.0;
    int const maxiter = 3;
    double const g = 2.0;

    ODEWithDirichletBCs ode;
    ode.known_solutions = {{{0}, {g}}};

    NumLib::BackwardEuler time_disc;
    DirichletNewtonODESystem ode_sys(process_id, ode, time_disc);

    auto linear_solver = createLinearSolver();
    RecordingResidualCriterion convergence_criterion;
    auto step_strategy = std::make_unique<RecordingStepStrategy>();
    auto const& recording_strategy = *step_strategy;
    NumLib::NonlinearSolver<NumLib::NonlinearSolverTag::Newton>
        nonlinear_solver(*linear_solver, maxiter, std::move(step_strategy));
    nonlinear_solver.setEquationSystem(ode_sys, convergence_criterion);

    GlobalVector x0(ode.getMatrixSpecifications(process_id).nrows);
    ODETraits<ODE1>::setIC(x0);  // x0[0] = 1 != g

    std::vector<GlobalVector*> xs{
        &NumLib::GlobalVectorProvider::provider.getVector(x0)};
    std::vector<GlobalVector*> xs_prev{
        &NumLib::GlobalVectorProvider::provider.getVector(x0)};

    time_disc.setInitialState(t0);
    time_disc.nextTimestep(t0 + dt, dt);

    nonlinear_solver.solve(xs, xs_prev, nullptr, process_id);

    auto const& deviations = recording_strategy.dirichlet_deviations;
    ASSERT_EQ(maxiter, static_cast<int>(deviations.size()));
    EXPECT_EQ(1.0, deviations[0]);  // |x0[0] - g|
    EXPECT_NEAR(0.0, deviations[1], 1e-14);
    EXPECT_NEAR(0.0, deviations[2], 1e-14);
    EXPECT_NEAR(g, (*xs.front())[0], 1e-14);

    ASSERT_EQ(maxiter,
              static_cast<int>(convergence_criterion.residuals.size()));
    // The BC violation drives the predictor step.
    EXPECT_NE(0.0, convergence_criterion.residuals[0][0]);
    EXPECT_NEAR(0.0, convergence_criterion.residuals[1][0], 1e-14);
    EXPECT_NEAR(0.0, convergence_criterion.residuals[2][0], 1e-14);

    NumLib::GlobalVectorProvider::provider.releaseVector(*xs.front());
    NumLib::GlobalVectorProvider::provider.releaseVector(*xs_prev.front());
}

/*! The merit at a trial point is the norm of the Dirichlet-free residual
 * P (r - J V) with V_k = x_trial,k - g_k, with the non-equilibrium initial
 * residuum subtracted.
 *
 * For ODE1 with backward Euler and dt = 1 the residual is
 * r = (x - x_prev) + K x with K = [[0, 1], [-1, 0]], and J = I + K. For
 * x_prev = (1, 0) and g = 2 at DOF 0:
 * - x_trial = (2, 1) fulfils the BC, V = 0, r = (2, -1), hence (0, -1);
 * - x_trial = (1.5, 1) violates it, V = (-0.5, 0), r = (1.5, -0.5),
 *   J V = (-0.5, 0.5), hence (0, -1) as well: the linearization is exact
 *   for this linear problem.
 * With r_neq = (5, 3) both become (0, -4).
 */
#ifndef USE_PETSC
TEST(NumLibDirichletResidual, EvaluateMeritAtTrialPoint)
#else
TEST(NumLibDirichletResidual, DISABLED_EvaluateMeritAtTrialPoint)
#endif
{
    int const process_id = 0;
    ODEWithDirichletBCs ode;
    ode.known_solutions = {{{0}, {2.0}}};

    NumLib::BackwardEuler time_disc;
    time_disc.setInitialState(0.0);
    time_disc.nextTimestep(1.0, 1.0);
    DirichletNewtonODESystem ode_sys(process_id, ode, time_disc);

    auto x = makeDirichletTestVector({2.0, 0.0});
    auto x_prev = makeDirichletTestVector({1.0, 0.0});
    std::vector<GlobalVector*> xs{&x};
    std::vector<GlobalVector*> xs_prev{&x_prev};
    ode_sys.computeKnownSolutions(x, process_id);

    auto x_trial = makeDirichletTestVector({2.0, 1.0});
    // A trial point violating the BC, e.g., after a damped step in the first
    // iteration of a time step.
    auto x_trial_violating = makeDirichletTestVector({1.5, 1.0});

    NumLib::NewtonStepContext const ctx{ode_sys, xs, xs_prev, process_id,
                                        nullptr};
    auto const merit = ctx.evaluateMerit(x_trial);
    ASSERT_TRUE(merit.has_value());
    EXPECT_NEAR(1.0, *merit, 1e-15);
    auto const merit_violating = ctx.evaluateMerit(x_trial_violating);
    ASSERT_TRUE(merit_violating.has_value());
    EXPECT_NEAR(1.0, *merit_violating, 1e-15);

    auto const r_neq = makeDirichletTestVector({5.0, 3.0});
    NumLib::NewtonStepContext const ctx_r_neq{ode_sys, xs, xs_prev, process_id,
                                              &r_neq};
    auto const merit_r_neq = ctx_r_neq.evaluateMerit(x_trial);
    ASSERT_TRUE(merit_r_neq.has_value());
    EXPECT_NEAR(4.0, *merit_r_neq, 1e-15);
    auto const merit_r_neq_violating =
        ctx_r_neq.evaluateMerit(x_trial_violating);
    ASSERT_TRUE(merit_r_neq_violating.has_value());
    EXPECT_NEAR(4.0, *merit_r_neq_violating, 1e-15);
}

//! A failing assembly at a trial point is reported, not thrown, s.t. a
//! globalization strategy can, e.g., backtrack.
#ifndef USE_PETSC
TEST(NumLibDirichletResidual, EvaluateMeritReportsAssemblyFailure)
#else
TEST(NumLibDirichletResidual, DISABLED_EvaluateMeritReportsAssemblyFailure)
#endif
{
    int const process_id = 0;
    ODEWithDirichletBCs ode;
    ode.assembly_fails = true;

    NumLib::BackwardEuler time_disc;
    time_disc.setInitialState(0.0);
    time_disc.nextTimestep(1.0, 1.0);
    DirichletNewtonODESystem ode_sys(process_id, ode, time_disc);

    auto x = makeDirichletTestVector({1.0, 0.0});
    std::vector<GlobalVector*> xs{&x};
    std::vector<GlobalVector*> xs_prev{&x};
    ode_sys.computeKnownSolutions(x, process_id);

    auto x_trial = makeDirichletTestVector({1.0, 1.0});

    NumLib::NewtonStepContext const ctx{ode_sys, xs, xs_prev, process_id,
                                        nullptr};
    EXPECT_NO_THROW({ EXPECT_FALSE(ctx.evaluateMerit(x_trial).has_value()); });
}

/*! If the set of Dirichlet DOFs changes between iterations of a solve, the
 * step strategy sees the deviation from the new Dirichlet value in that
 * iteration, and the full step imposes it.
 */
#ifndef USE_PETSC
TEST(NumLibDirichletResidual, DeviationAfterChangedSetOfDirichletDofs)
#else
TEST(NumLibDirichletResidual, DISABLED_DeviationAfterChangedSetOfDirichletDofs)
#endif
{
    int const process_id = 0;
    int const maxiter = 4;

    ODEWithDirichletBCs ode;
    // Iterations 1, 2: DOF 0 is fixed, iterations 3, 4: DOF 1 is fixed.
    ode.known_solutions_per_call = {
        {{{0}, {2.0}}}, {{{0}, {2.0}}}, {{{1}, {0.5}}}, {{{1}, {0.5}}}};

    NumLib::BackwardEuler time_disc;
    DirichletNewtonODESystem ode_sys(process_id, ode, time_disc);

    auto linear_solver = createLinearSolver();
    RecordingResidualCriterion convergence_criterion;
    auto step_strategy = std::make_unique<RecordingStepStrategy>();
    auto const& recording_strategy = *step_strategy;
    NumLib::NonlinearSolver<NumLib::NonlinearSolverTag::Newton>
        nonlinear_solver(*linear_solver, maxiter, std::move(step_strategy));
    nonlinear_solver.setEquationSystem(ode_sys, convergence_criterion);

    GlobalVector x0(ode.getMatrixSpecifications(process_id).nrows);
    ODETraits<ODE1>::setIC(x0);

    std::vector<GlobalVector*> xs{
        &NumLib::GlobalVectorProvider::provider.getVector(x0)};
    std::vector<GlobalVector*> xs_prev{
        &NumLib::GlobalVectorProvider::provider.getVector(x0)};

    time_disc.setInitialState(0.0);
    time_disc.nextTimestep(1.0, 1.0);

    nonlinear_solver.solve(xs, xs_prev, nullptr, process_id);

    // The first iteration imposes DOF 0: x0 = (1, 0) -> (2, 2). The third
    // one imposes DOF 1.
    auto const& deviations = recording_strategy.dirichlet_deviations;
    ASSERT_EQ(maxiter, static_cast<int>(deviations.size()));
    EXPECT_EQ(1.0, deviations[0]);
    EXPECT_NEAR(0.0, deviations[1], 1e-14);
    EXPECT_NEAR(1.5, deviations[2], 1e-14);
    EXPECT_NEAR(0.0, deviations[3], 1e-14);

    // After each predictor the fixed DOF's residual vanishes.
    ASSERT_EQ(maxiter,
              static_cast<int>(convergence_criterion.residuals.size()));
    EXPECT_NEAR(0.0, convergence_criterion.residuals[1][0], 1e-14);
    EXPECT_NEAR(0.0, convergence_criterion.residuals[3][1], 1e-14);
    EXPECT_NEAR(0.5, (*xs.front())[1], 1e-14);

    NumLib::GlobalVectorProvider::provider.releaseVector(*xs.front());
    NumLib::GlobalVectorProvider::provider.releaseVector(*xs_prev.front());
}

/*! With dirichlet_free_residual_check, the convergence criterion sees
 * P (res - J V), V_k = x_k - g_k: the residual linearized at the state
 * fulfilling the BCs, without the Dirichlet entries.
 *
 * ODE1 with backward Euler, dt = 1, x0 = x_prev = (1, 0), g = 2 at DOF 0:
 * res = (x - x_prev) + K x = (0, -1), J = I + K = [[1, 1], [-1, 1]],
 * V = (-1, 0), hence res - J V = (1, -2) and the checked vector is (0, -2).
 * The default check sees Eigen's residual after applying the BCs instead,
 * (c v, res_1 - J_10 v) = (-1, -2) with c = J_00 = 1.
 * After the full step x = (2, 2) fulfils the BC, and the checked vector is
 * (0, 0) since the problem is linear.
 * The merit passed to the step strategy is the norm of the Dirichlet-free
 * residual, (2, 0), regardless of the option.
 */
#ifndef USE_PETSC
TEST(NumLibDirichletResidual, DirichletFreeResidualCheck)
#else
TEST(NumLibDirichletResidual, DISABLED_DirichletFreeResidualCheck)
#endif
{
    int const process_id = 0;
    int const maxiter = 2;

    std::vector<std::vector<double>> merits;
    auto checked_residuals = [&](bool const dirichlet_free)
    {
        ODEWithDirichletBCs ode;
        ode.known_solutions = {{{0}, {2.0}}};

        NumLib::BackwardEuler time_disc;
        DirichletNewtonODESystem ode_sys(process_id, ode, time_disc);

        auto linear_solver = createLinearSolver();
        RecordingResidualCriterion convergence_criterion;
        auto step_strategy = std::make_unique<RecordingStepStrategy>();
        auto const& recording_strategy = *step_strategy;
        NumLib::NonlinearSolver<NumLib::NonlinearSolverTag::Newton>
            nonlinear_solver(*linear_solver, maxiter, std::move(step_strategy));
        nonlinear_solver.setDirichletFreeResidualCheck(dirichlet_free);
        nonlinear_solver.setEquationSystem(ode_sys, convergence_criterion);

        GlobalVector x0(ode.getMatrixSpecifications(process_id).nrows);
        ODETraits<ODE1>::setIC(x0);  // (1, 0)

        std::vector<GlobalVector*> xs{
            &NumLib::GlobalVectorProvider::provider.getVector(x0)};
        std::vector<GlobalVector*> xs_prev{
            &NumLib::GlobalVectorProvider::provider.getVector(x0)};

        time_disc.setInitialState(0.0);
        time_disc.nextTimestep(1.0, 1.0);

        nonlinear_solver.solve(xs, xs_prev, nullptr, process_id);
        merits.push_back(recording_strategy.merits);

        NumLib::GlobalVectorProvider::provider.releaseVector(*xs.front());
        NumLib::GlobalVectorProvider::provider.releaseVector(*xs_prev.front());
        return std::move(convergence_criterion.residuals);
    };

    auto const dirichlet_free = checked_residuals(true);
    ASSERT_EQ(maxiter, static_cast<int>(dirichlet_free.size()));
    EXPECT_EQ(0.0, dirichlet_free[0][0]);
    EXPECT_NEAR(-2.0, dirichlet_free[0][1], 1e-14);
    EXPECT_EQ(0.0, dirichlet_free[1][0]);
    EXPECT_NEAR(0.0, dirichlet_free[1][1], 1e-14);

    auto const with_dirichlet_rows = checked_residuals(false);
    ASSERT_EQ(maxiter, static_cast<int>(with_dirichlet_rows.size()));
    EXPECT_NEAR(-1.0, with_dirichlet_rows[0][0], 1e-14);
    EXPECT_NEAR(-2.0, with_dirichlet_rows[0][1], 1e-14);

    for (auto const& run_merits : merits)
    {
        ASSERT_EQ(maxiter, static_cast<int>(run_merits.size()));
        EXPECT_NEAR(2.0, run_merits[0], 1e-14);
        EXPECT_NEAR(0.0, run_merits[1], 1e-14);
    }
}

#ifdef USE_PETSC
/*! For PETSc SNES the residual at a Dirichlet DOF is the deviation of the
 * iterate from the BC value, s.t. SNES notices if one of its steps, e.g., a
 * trust region step, moves a Dirichlet DOF. The iterate itself must not be
 * modified, other residual entries neither.
 *
 * Unlike the other tests here this one runs with PETSc only, because the
 * SNES path is not used otherwise. It provides the sparsity pattern needed
 * for PETSc matrices (cf. issue #1989).
 */
TEST(NumLibDirichletResidual, PETScSNESResidualIsDeviationFromBC)
{
    int const process_id = 0;
    double const g = 2.0;

    auto const sparsity_pattern =
        PETScSparsityPatternUtils::denseBlockSparsityPattern(2, 2, 2, 0);
    ODEWithDirichletBCs ode;
    ode.sparsity_pattern = &sparsity_pattern;
    ode.known_solutions = {{{0}, {g}}};

    NumLib::BackwardEuler time_disc;
    DirichletNewtonODESystem ode_sys(process_id, ode, time_disc);

    auto J = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(
        ode.getMatrixSpecifications(process_id));
    MathLib::setMatrix(*J, {4.0, 1.0, 1.0, 3.0});
    MathLib::LinAlg::finalizeAssembly(*J);

    // x[0] deviates from the BC value g by -0.5.
    auto x = makeDirichletTestVector({1.5, 0.25});
    auto res = makeDirichletTestVector({7.0, 8.0});
    auto scratch = makeDirichletTestVector({0.0, 0.0});

    ode_sys.computeKnownSolutions(x, process_id);
    MathLib::LinAlg::setLocalAccessibleVector(x);
    ode_sys.applyKnownSolutionsPETScSNES(*J, res, x, scratch);

    MathLib::LinAlg::setLocalAccessibleVector(res);
    EXPECT_EQ(1.5 - g, res[0]);
    EXPECT_EQ(8.0, res[1]);

    MathLib::LinAlg::setLocalAccessibleVector(x);
    EXPECT_EQ(1.5, x[0]);
    EXPECT_EQ(0.25, x[1]);
}
#endif
