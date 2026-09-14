// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>
#include <optional>

#include "CreateTestLinearSolver.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/Exceptions.h"
#include "NumLib/NumericsConfig.h"
#include "NumLib/ODESolver/ConvergenceCriterionDeltaX.h"
#include "NumLib/ODESolver/FixedDampingStrategy.h"
#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "ODEs.h"

namespace
{
/*! ODE1 whose assemblies throw NumLib::AssemblyException once, if so
 * requested.
 *
 * It models a local assembler that cannot compute its contribution, e.g.,
 * because a local nonlinear solver did not converge. Later assemblies
 * succeed, s.t. a repeated time step can converge.
 *
 * ODE1 is tagged Newton, s.t. it can be used with both nonlinear solvers.
 */
class ODEAssemblyFails final : public ODE1
{
public:
    //! \param failing_assembly One-based number of the assembly that throws,
    //! std::nullopt if no assembly throws.
    explicit ODEAssemblyFails(std::optional<int> const failing_assembly)
        : failing_assembly_(failing_assembly)
    {
    }

    void assemble(const double t, double const dt,
                  std::vector<GlobalVector*> const& x,
                  std::vector<GlobalVector*> const& x_prev,
                  int const process_id, GlobalMatrix& M, GlobalMatrix& K,
                  GlobalVector& b) override
    {
        throwIfRequested();
        ODE1::assemble(t, dt, x, x_prev, process_id, M, K, b);
    }

    void assembleWithJacobian(const double t, double const dt,
                              std::vector<GlobalVector*> const& x_curr,
                              std::vector<GlobalVector*> const& x_prev,
                              int const process_id, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        throwIfRequested();
        ODE1::assembleWithJacobian(t, dt, x_curr, x_prev, process_id, b, Jac);
    }

    /*! ODE1 is linear, and both solvers stop after a single iteration on a
     * linear system. Reporting the system as nonlinear makes them iterate
     * until the convergence criterion is met, s.t. an assembly after the
     * first one is reached at all.
     */
    bool isLinear() const override { return false; }

    //! Number of times an assembly has been entered.
    int number_of_assemblies = 0;

private:
    void throwIfRequested()
    {
        ++number_of_assemblies;

        if (number_of_assemblies == failing_assembly_)
        {
            throw NumLib::AssemblyException("ODEAssemblyFails: by design.");
        }
    }

    std::optional<int> const failing_assembly_;
};

/*! Runs a single nonlinear solve of \c ode over one time step.
 *
 * If the assembly of \c ode throws, this exercises the AssemblyException
 * handling of the nonlinear solver.
 */
template <NumLib::NonlinearSolverTag NLTag>
NumLib::NonlinearSolverStatus solveOneStep(ODEAssemblyFails& ode)
{
    int const process_id = 0;
    double const t0 = 0.0;
    double const dt = 1.0;
    int const maxiter = 20;

    NumLib::BackwardEuler time_disc;
    NumLib::TimeDiscretizedODESystem<ODEAssemblyFails::ODETag, NLTag> ode_sys(
        process_id, ode, time_disc);

    auto linear_solver = createLinearSolver();
    auto convergence_criterion =
        std::make_unique<NumLib::ConvergenceCriterionDeltaX>(
            1e-9, std::nullopt, MathLib::VecNormType::NORM2);
    auto nonlinear_solver = [&linear_solver, maxiter]
    {
        using NLSolver = NumLib::NonlinearSolver<NLTag>;
        if constexpr (NLTag == NumLib::NonlinearSolverTag::Newton)
        {
            return std::make_unique<NLSolver>(
                *linear_solver, maxiter,
                std::make_unique<NumLib::FixedDampingStrategy>(1.0));
        }
        else
        {
            int const anderson_depth = 0;
            return std::make_unique<NLSolver>(*linear_solver, maxiter,
                                              anderson_depth, 1.0);
        }
    }();
    nonlinear_solver->setEquationSystem(ode_sys, *convergence_criterion);

    GlobalVector x0(ode.getMatrixSpecifications(process_id).nrows);
    ODETraits<ODE1>::setIC(x0);

    std::vector<GlobalVector*> xs{
        &NumLib::GlobalVectorProvider::provider.getVector(x0)};
    std::vector<GlobalVector*> xs_prev{
        &NumLib::GlobalVectorProvider::provider.getVector(x0)};

    time_disc.setInitialState(t0);
    MathLib::LinAlg::copy(*xs.front(), *xs_prev.front());
    time_disc.nextTimestep(t0 + dt, dt);

    auto const status =
        nonlinear_solver->solve(xs, xs_prev, nullptr, process_id);

    NumLib::GlobalVectorProvider::provider.releaseVector(*xs.front());
    NumLib::GlobalVectorProvider::provider.releaseVector(*xs_prev.front());

    return status;
}

/*! Checks that a failing assembly aborts the nonlinear iteration.
 *
 * The AssemblyException must not escape the nonlinear solver, but be reported
 * as a non-converged time step, s.t. the time stepping algorithm can repeat
 * the step with a smaller time step size.
 *
 * \c failing_assembly selects whether the solve is abandoned right away or
 * only after an iteration has already moved the iterate: the reported number
 * of iterations is the one the assembly failed in either way.
 */
template <NumLib::NonlinearSolverTag NLTag>
void expectAssemblyExceptionAbortsIteration(int const failing_assembly)
{
    ODEAssemblyFails ode{failing_assembly};

    NumLib::NonlinearSolverStatus status;
    ASSERT_NO_THROW(status = solveOneStep<NLTag>(ode));

    EXPECT_FALSE(status.error_norms_met);
    // The iteration the assembly actually failed in is reported, not _maxiter.
    EXPECT_EQ(failing_assembly, status.number_iterations);
    EXPECT_EQ(failing_assembly, ode.number_of_assemblies);
}

//! Checks that the very same setup converges without a failing assembly.
//! Guards against the exception handling breaking the regular path.
template <NumLib::NonlinearSolverTag NLTag>
void expectConvergenceWithoutAssemblyException()
{
    ODEAssemblyFails ode{std::nullopt};

    NumLib::NonlinearSolverStatus status;
    ASSERT_NO_THROW(status = solveOneStep<NLTag>(ode));

    EXPECT_TRUE(status.error_norms_met);
}
}  // namespace

//! The Picard solver must survive an AssemblyException.
// Disabled for PETSc together with the other ODE system tests, see issue
// #1989: the ODE systems used for testing do not provide a sparsity pattern,
// which PETSc requires for matrix preallocation.
#ifndef USE_PETSC
TEST(NumLibNonlinearSolver, PicardAssemblyExceptionAbortsIteration)
#else
TEST(NumLibNonlinearSolver, DISABLED_PicardAssemblyExceptionAbortsIteration)
#endif
{
    expectAssemblyExceptionAbortsIteration<NumLib::NonlinearSolverTag::Picard>(
        1);
}

#ifndef USE_PETSC
TEST(NumLibNonlinearSolver, PicardWithoutAssemblyExceptionConverges)
#else
TEST(NumLibNonlinearSolver, DISABLED_PicardWithoutAssemblyExceptionConverges)
#endif
{
    expectConvergenceWithoutAssemblyException<
        NumLib::NonlinearSolverTag::Picard>();
}

//! Same as above for the Newton solver, which has caught the exception since
//! it was introduced.
#ifndef USE_PETSC
TEST(NumLibNonlinearSolver, NewtonAssemblyExceptionAbortsIteration)
#else
TEST(NumLibNonlinearSolver, DISABLED_NewtonAssemblyExceptionAbortsIteration)
#endif
{
    expectAssemblyExceptionAbortsIteration<NumLib::NonlinearSolverTag::Newton>(
        1);
}

#ifndef USE_PETSC
TEST(NumLibNonlinearSolver, NewtonWithoutAssemblyExceptionConverges)
#else
TEST(NumLibNonlinearSolver, DISABLED_NewtonWithoutAssemblyExceptionConverges)
#endif
{
    expectConvergenceWithoutAssemblyException<
        NumLib::NonlinearSolverTag::Newton>();
}

//! The abort must also work once an iteration has already moved the iterate,
//! and still report the iteration the assembly failed in.
#ifndef USE_PETSC
TEST(NumLibNonlinearSolver, PicardAssemblyExceptionAbortsIterationAfterProgress)
#else
TEST(NumLibNonlinearSolver,
     DISABLED_PicardAssemblyExceptionAbortsIterationAfterProgress)
#endif
{
    expectAssemblyExceptionAbortsIteration<NumLib::NonlinearSolverTag::Picard>(
        2);
}

#ifndef USE_PETSC
TEST(NumLibNonlinearSolver, NewtonAssemblyExceptionAbortsIterationAfterProgress)
#else
TEST(NumLibNonlinearSolver,
     DISABLED_NewtonAssemblyExceptionAbortsIterationAfterProgress)
#endif
{
    expectAssemblyExceptionAbortsIteration<NumLib::NonlinearSolverTag::Newton>(
        2);
}
