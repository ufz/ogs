// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/ODESolver/ConvergenceCriterion.h"
#include "NumLib/ODESolver/NonlinearSolverStatus.h"
#include "NumLib/StaggeredCoupling/StaggeredCoupling.h"

namespace
{
using CouplingNodeVariant = NumLib::RootCouplingNode::CouplingNodeVariant;

/// Convergence criterion whose delta-x check never passes. It drives the
/// coupling loop into exhausting its iteration budget.
class NeverSatisfiedCriterion final : public NumLib::ConvergenceCriterion
{
public:
    NeverSatisfiedCriterion()
        : NumLib::ConvergenceCriterion(MathLib::VecNormType::NORM2)
    {
    }

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return false; }

    void checkDeltaX(GlobalVector const& /*minus_delta_x*/,
                     GlobalVector const& /*x*/) override
    {
        _satisfied = false;
    }

    void checkResidual(GlobalVector const& /*residual*/) override {}
};

/// Counterpart of NeverSatisfiedCriterion: the delta-x check always passes.
class AlwaysSatisfiedCriterion final : public NumLib::ConvergenceCriterion
{
public:
    AlwaysSatisfiedCriterion()
        : NumLib::ConvergenceCriterion(MathLib::VecNormType::NORM2)
    {
    }

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return false; }

    void checkDeltaX(GlobalVector const& /*minus_delta_x*/,
                     GlobalVector const& /*x*/) override
    {
    }

    void checkResidual(GlobalVector const& /*residual*/) override {}
};

struct TestProcessData
{
    NumLib::NonlinearSolverStatus nonlinear_solver_status = {true, 1};
};

struct TestOutput
{
};

template <typename Criterion>
std::vector<CouplingNodeVariant> createCouplingNodes()
{
    std::vector<CouplingNodeVariant> coupling_nodes;
    for (int process_id = 0; process_id < 2; ++process_id)
    {
        coupling_nodes.emplace_back(
            NumLib::CouplingNode{"process_" + std::to_string(process_id),
                                 std::make_unique<Criterion>(), 1, process_id});
    }
    return coupling_nodes;
}

std::vector<std::unique_ptr<TestProcessData>> createProcessData()
{
    std::vector<std::unique_ptr<TestProcessData>> per_process_data;
    per_process_data.emplace_back(std::make_unique<TestProcessData>());
    per_process_data.emplace_back(std::make_unique<TestProcessData>());
    return per_process_data;
}

/// Solves nothing; it only reports the given status and changes the solution
/// so that the coupling convergence check operates on a non-zero increment.
auto createProcessSolver(bool const nonlinear_solver_converged)
{
    return [nonlinear_solver_converged, value = 0.0](
               std::vector<GlobalVector*>& xs,
               std::vector<GlobalVector*> const& /*xs_prev*/,
               std::size_t const /*timestep*/, double const /*t*/,
               double const /*dt*/, TestProcessData const& /*process_data*/,
               std::vector<TestOutput> const& /*outputs*/) mutable
    {
        value += 1.0;
        for (auto* x : xs)
        {
            MathLib::LinAlg::set(*x, value);
        }
        return NumLib::NonlinearSolverStatus{nonlinear_solver_converged, 1};
    };
}

/// Stands in for the production TimeLoop: it owns the solution vectors, the
/// per-process data and the outputs, and drives one staggered coupling step
/// (mirrors TimeLoop::solveCoupledEquationSystemsByStaggeredScheme).
struct TimeLoopMockup
{
    TimeLoopMockup()
    {
        for (auto* x : {&x0, &x1, &x0_prev, &x1_prev})
        {
            MathLib::LinAlg::set(*x, 0.0);
        }
    }

    NumLib::NonlinearSolverStatus solveCoupledTimeStep(
        NumLib::StaggeredCoupling& coupling,
        bool const nonlinear_solver_converged)
    {
        coupling.initializeCoupledSolutions(xs);
        return coupling.execute<TestProcessData, TestOutput>(
            /*t=*/1.0, /*dt=*/1.0, /*timestep_id=*/1, xs, xs_prev,
            per_process_data, outputs,
            createProcessSolver(nonlinear_solver_converged));
    }

    GlobalVector x0{3}, x1{3}, x0_prev{3}, x1_prev{3};
    std::vector<GlobalVector*> xs = {&x0, &x1};
    std::vector<GlobalVector*> const xs_prev = {&x0_prev, &x1_prev};
    std::vector<std::unique_ptr<TestProcessData>> per_process_data =
        createProcessData();
    std::vector<TestOutput> const outputs;
};
}  // namespace

class NumLibStaggeredCoupling : public ::testing::Test
{
protected:
    TimeLoopMockup time_loop;
};

#ifndef USE_PETSC
TEST_F(NumLibStaggeredCoupling, RejectsTimeStepIfCouplingDoesNotConverge)
#else
TEST_F(NumLibStaggeredCoupling,
       DISABLED_RejectsTimeStepIfCouplingDoesNotConverge)
#endif
{
    int const max_coupling_iterations = 3;
    NumLib::StaggeredCoupling coupling{
        max_coupling_iterations,
        createCouplingNodes<NeverSatisfiedCriterion>()};

    auto const status = time_loop.solveCoupledTimeStep(
        coupling, /*nonlinear_solver_converged=*/true);

    EXPECT_EQ(max_coupling_iterations,
              coupling.lastNumberOfCouplingIterations());

    // The per-process nonlinear solvers converged on their own, but the
    // coupling did not. The time loop must not accept such a solution.
    EXPECT_FALSE(status.error_norms_met);
    for (auto const& process_data : time_loop.per_process_data)
    {
        EXPECT_FALSE(process_data->nonlinear_solver_status.error_norms_met);
    }
}

#ifndef USE_PETSC
TEST_F(NumLibStaggeredCoupling, AcceptsTimeStepIfCouplingConverges)
#else
TEST_F(NumLibStaggeredCoupling, DISABLED_AcceptsTimeStepIfCouplingConverges)
#endif
{
    NumLib::StaggeredCoupling coupling{
        /*global_coupling_max_iterations=*/3,
        createCouplingNodes<AlwaysSatisfiedCriterion>()};

    auto const status = time_loop.solveCoupledTimeStep(
        coupling, /*nonlinear_solver_converged=*/true);

    EXPECT_TRUE(status.error_norms_met);
    for (auto const& process_data : time_loop.per_process_data)
    {
        EXPECT_TRUE(process_data->nonlinear_solver_status.error_norms_met);
    }
}

#ifndef USE_PETSC
TEST_F(NumLibStaggeredCoupling, RejectsTimeStepIfNonlinearSolverFails)
#else
TEST_F(NumLibStaggeredCoupling, DISABLED_RejectsTimeStepIfNonlinearSolverFails)
#endif
{
    NumLib::StaggeredCoupling coupling{
        /*global_coupling_max_iterations=*/3,
        createCouplingNodes<AlwaysSatisfiedCriterion>()};

    auto const status = time_loop.solveCoupledTimeStep(
        coupling, /*nonlinear_solver_converged=*/false);

    // A failed nonlinear solver leaves the coupling loop immediately.
    EXPECT_EQ(1, coupling.lastNumberOfCouplingIterations());
    EXPECT_FALSE(status.error_norms_met);
}
