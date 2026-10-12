// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include "NumLib/StaggeredCoupling/StaggeredCoupling.h"
#include "Tests/NumLib/StaggeredCouplingTestDoubles.h"

namespace
{

class NumLibStaggeredCoupling : public ::testing::Test
{
protected:
    StaggeredCouplingTestDoubles::CoupledSystem system;
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
        StaggeredCouplingTestDoubles::createCouplingNodes<
            StaggeredCouplingTestDoubles::NeverSatisfiedCriterion>()};

    auto const status =
        system.solveOneStep(coupling, /*nonlinear_solver_converged=*/true);

    EXPECT_EQ(max_coupling_iterations,
              coupling.lastNumberOfCouplingIterations());

    // The per-process nonlinear solvers converged on their own, but the
    // coupling did not. The time loop must not accept such a solution.
    EXPECT_FALSE(status.error_norms_met);
    for (auto const& process_data : system.per_process_data)
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
        StaggeredCouplingTestDoubles::createCouplingNodes<
            StaggeredCouplingTestDoubles::AlwaysSatisfiedCriterion>()};

    auto const status =
        system.solveOneStep(coupling, /*nonlinear_solver_converged=*/true);

    EXPECT_TRUE(status.error_norms_met);
    for (auto const& process_data : system.per_process_data)
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
        StaggeredCouplingTestDoubles::createCouplingNodes<
            StaggeredCouplingTestDoubles::AlwaysSatisfiedCriterion>()};

    auto const status =
        system.solveOneStep(coupling, /*nonlinear_solver_converged=*/false);

    // A failed nonlinear solver leaves the coupling loop immediately.
    EXPECT_EQ(1, coupling.lastNumberOfCouplingIterations());
    EXPECT_FALSE(status.error_norms_met);
}
}  // namespace
