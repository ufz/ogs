// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#if defined(USE_LIS)
#include "MathLib/LinAlg/EigenLis/LinearSolverOptionsParser.h"
#elif defined(USE_PETSC)
#include "MathLib/LinAlg/PETSc/LinearSolverOptionsParser.h"
#else
#include "MathLib/LinAlg/Eigen/LinearSolverOptionsParser.h"
#endif
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/Exceptions.h"
#include "NumLib/NumericsConfig.h"
#include "NumLib/ODESolver/ConvergenceCriterionDeltaX.h"
#include "NumLib/ODESolver/FixedDampingStrategy.h"
#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"

namespace
{
// An ODE system whose assembly reports a state it cannot assemble, the way the
// drift-flux closure of the WellboreSimulator process and the MFront and
// Coulomb material models do. Tagged for Newton so that it satisfies the
// interface of both nonlinear solvers.
//
// The number of assemblies that succeed before the bad state is reported is a
// parameter, because the two cases exercise different things. Aborting on the
// first assembly is the degenerate case, where no iterate has been computed
// yet. Aborting after a few successful ones is what the closures actually do:
// they only reach the inadmissible region once the iterate has moved into it,
// so the solver has to abandon a partially progressed solution rather than an
// untouched one.
class AbortingODE final
    : public NumLib::ODESystem<
          NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
          NumLib::NonlinearSolverTag::Newton>
{
private:
    // The system is dense, and PETSc preallocates from the pattern rather than
    // from the matrix, so it has to be given even though this assembly never
    // writes an entry: it throws before it does.
    static GlobalSparsityPattern densePattern()
    {
        GlobalSparsityPattern pattern;
#ifdef USE_PETSC
        pattern.row_ptr.push_back(0);
        for (std::size_t row = 0; row < N; ++row)
        {
            for (std::size_t column = 0; column < N; ++column)
            {
                pattern.col_idx.push_back(static_cast<PetscInt>(column));
            }
            pattern.row_ptr.push_back(pattern.row_ptr.back() +
                                      static_cast<PetscInt>(N));
        }
#else
        pattern.number_non_zeros_per_row.assign(
            N, static_cast<GlobalIndexType>(N));
#endif
        return pattern;
    }

    GlobalSparsityPattern const sparsity_pattern_ = densePattern();

    int const assemblies_before_abort_;
    int assemblies_ = 0;

    // Counts this assembly and reports the bad state once the configured
    // number of successful ones has been spent.
    void abortOnceTheIterateHasMoved()
    {
        if (assemblies_++ >= assemblies_before_abort_)
        {
            throw NumLib::AssemblyException(
                "Test assembly reports a bad state.");
        }
    }

    // Grows by one per assembly, so that no two consecutive iterates coincide
    // and the solve cannot converge before it aborts.
    double movingRightHandSide() const
    {
        return static_cast<double>(assemblies_);
    }

public:
    explicit AbortingODE(int const assemblies_before_abort = 0)
        : assemblies_before_abort_(assemblies_before_abort)
    {
    }

    void preAssemble(double const /*t*/, double const /*dt*/,
                     GlobalVector const& /*x*/) override
    {
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<GlobalVector*> const& /*x*/,
                  std::vector<GlobalVector*> const& /*x_prev*/,
                  int const /*process_id*/, GlobalMatrix& M, GlobalMatrix& K,
                  GlobalVector& b) override
    {
        abortOnceTheIterateHasMoved();

        // The identity in K with an empty M leaves the Picard matrix as the
        // identity, so the iterate is whatever the right-hand side says. The
        // right-hand side moves by one in every component per assembly, which
        // keeps the increment away from the convergence criterion until the
        // abort arrives.
        MathLib::setMatrix(M, {0.0, 0.0, 0.0, 0.0});
        MathLib::setMatrix(K, {1.0, 0.0, 0.0, 1.0});
        MathLib::setVector(b, {movingRightHandSide(), movingRightHandSide()});
    }

    void assembleWithJacobian(double const /*t*/, double const /*dt*/,
                              std::vector<GlobalVector*> const& /*x*/,
                              std::vector<GlobalVector*> const& /*x_prev*/,
                              int const /*process_id*/, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        abortOnceTheIterateHasMoved();

        // As above: an identity Jacobian and a residual that grows by one per
        // assembly, so every Newton step moves the iterate by at least one.
        MathLib::setMatrix(Jac, {1.0, 0.0, 0.0, 1.0});
        MathLib::setVector(b, {movingRightHandSide(), movingRightHandSide()});
    }

    MathLib::MatrixSpecifications getMatrixSpecifications(
        int const /*process_id*/) const override
    {
        return {N, N, nullptr, &sparsity_pattern_};
    }

    bool isLinear() const override { return false; }

    bool requiresNormalization() const override { return false; }

    static constexpr std::size_t N = 2;
};

std::unique_ptr<GlobalLinearSolver> createLinearSolver()
{
#if defined(USE_PETSC)
    return std::make_unique<GlobalLinearSolver>(
        "", "-ksp_type bcgs -pc_type sor -ksp_rtol 1e-24 -ksp_max_it 100");
#else
    auto const solver_options =
        MathLib::LinearSolverOptionsParser<GlobalLinearSolver>{}
            .parseNameAndOptions("", nullptr);
    return std::make_unique<GlobalLinearSolver>(std::get<0>(solver_options),
                                                std::get<1>(solver_options));
#endif
}

// Iteration limit of the solver under test. The retreat reports it as the
// number of iterations spent, which is what the time stepper reads.
constexpr unsigned maxiter = 20;

// Runs a single nonlinear solve on the aborting ODE system and returns what
// the solver reports for it. The assembly aborts after \c
// assemblies_before_abort successful assemblies.
template <NumLib::NonlinearSolverTag NLTag>
NumLib::NonlinearSolverStatus solveOnce(int const assemblies_before_abort = 0)
{
    int const process_id = 0;

    AbortingODE ode{assemblies_before_abort};
    NumLib::BackwardEuler time_disc;
    NumLib::TimeDiscretizedODESystem<AbortingODE::ODETag, NLTag> ode_sys(
        process_id, ode, time_disc);

    auto linear_solver = createLinearSolver();
    auto convergence_criterion =
        std::make_unique<NumLib::ConvergenceCriterionDeltaX>(
            1e-9, std::nullopt, MathLib::VecNormType::NORM2);

    std::unique_ptr<NumLib::NonlinearSolver<NLTag>> nonlinear_solver;
    if constexpr (NLTag == NumLib::NonlinearSolverTag::Newton)
    {
        nonlinear_solver = std::make_unique<NumLib::NonlinearSolver<NLTag>>(
            *linear_solver, maxiter,
            std::make_unique<NumLib::FixedDampingStrategy>(1.0));
    }
    else
    {
        nonlinear_solver = std::make_unique<NumLib::NonlinearSolver<NLTag>>(
            *linear_solver, maxiter, 1.0);
    }
    nonlinear_solver->setEquationSystem(ode_sys, *convergence_criterion);

    GlobalVector x0(AbortingODE::N);
    MathLib::setVector(x0, {0.0, 0.0});
    MathLib::LinAlg::finalizeAssembly(x0);

    std::vector<GlobalVector*> xs{
        &NumLib::GlobalVectorProvider::provider.getVector(x0)};
    std::vector<GlobalVector*> xs_prev{
        &NumLib::GlobalVectorProvider::provider.getVector(x0)};

    time_disc.setInitialState(0.);
    time_disc.nextTimestep(1., 1.);

    auto const status =
        nonlinear_solver->solve(xs, xs_prev, nullptr, process_id);

    NumLib::GlobalVectorProvider::provider.releaseVector(*xs[0]);
    NumLib::GlobalVectorProvider::provider.releaseVector(*xs_prev[0]);

    return status;
}
}  // namespace

// An AssemblyException asks for the nonlinear iteration to be abandoned and
// the time step to be repeated, so it must not leave the nonlinear solver. The
// solver reports an unconverged iteration instead, which the time loop answers
// with a smaller time step.
//
// The reported iteration count is part of that answer, not a detail: the
// iteration number based time stepper interpolates the next step size from it,
// so an abort has to look like an exhausted iteration rather than like one
// cheap iteration that failed.
TEST(NumLibNonlinearSolver, NewtonRetreatsOnAbortedAssembly)
{
    NumLib::NonlinearSolverStatus status;

    ASSERT_NO_THROW(status = solveOnce<NumLib::NonlinearSolverTag::Newton>());
    EXPECT_FALSE(status.error_norms_met);
    EXPECT_EQ(static_cast<int>(maxiter), status.number_iterations);
}

TEST(NumLibNonlinearSolver, PicardRetreatsOnAbortedAssembly)
{
    NumLib::NonlinearSolverStatus status;

    ASSERT_NO_THROW(status = solveOnce<NumLib::NonlinearSolverTag::Picard>());
    EXPECT_FALSE(status.error_norms_met);
    EXPECT_EQ(static_cast<int>(maxiter), status.number_iterations);
}

// The realistic case: the closures reach the inadmissible region only after
// the iterate has moved into it, so the abort has to abandon a partially
// progressed solve. The reported status must not depend on how far it got,
// or an abort late in the iteration would look to the time stepper like a
// solve that nearly converged and would earn a step size it cannot handle.
TEST(NumLibNonlinearSolver, NewtonRetreatsOnAssemblyAbortedAfterProgress)
{
    NumLib::NonlinearSolverStatus status;

    ASSERT_NO_THROW(status = solveOnce<NumLib::NonlinearSolverTag::Newton>(3));
    EXPECT_FALSE(status.error_norms_met);
    EXPECT_EQ(static_cast<int>(maxiter), status.number_iterations);
}

TEST(NumLibNonlinearSolver, PicardRetreatsOnAssemblyAbortedAfterProgress)
{
    NumLib::NonlinearSolverStatus status;

    ASSERT_NO_THROW(status = solveOnce<NumLib::NonlinearSolverTag::Picard>(3));
    EXPECT_FALSE(status.error_norms_met);
    EXPECT_EQ(static_cast<int>(maxiter), status.number_iterations);
}
