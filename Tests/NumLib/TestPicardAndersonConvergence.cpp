// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>

#include "CreateTestLinearSolver.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/NumericsConfig.h"
#include "NumLib/ODESolver/ConvergenceCriterionDeltaX.h"
#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "ODEs.h"

namespace
{
/*! The Picard map g(x) = G x + c with G = [[1, -1/2], [1, 0]], c = (1, 0),
 * and the fixed point x* = (2, 2).
 *
 * With M = 0 and K = I the Picard step solves x_new = b(x) = g(x).
 *
 * From x_0 = 0 the first two steps are f_0 = (1, 0) and f_1 = (1, 1). Since
 * <f_0, f_1> = |f_0|^2, Anderson acceleration of depth 2 picks the weights
 * theta = (1, 0), whose mixed iterate x_0 + f_0 = x_1 coincides with the
 * iterate x_1 = (1, 0), although x_1 is not a fixed point.
 */
class ODEMixingReproducesIterate final : public ODE1
{
public:
    void assemble(const double /*t*/, double const /*dt*/,
                  std::vector<GlobalVector*> const& x,
                  std::vector<GlobalVector*> const& /*x_prev*/,
                  int const process_id, GlobalMatrix& M, GlobalMatrix& K,
                  GlobalVector& b) override
    {
        auto const& x_curr = *x[process_id];
        MathLib::LinAlg::setLocalAccessibleVector(x_curr);
        double const x0 = x_curr.get(0);
        double const x1 = x_curr.get(1);

        MathLib::setMatrix(M, {0.0, 0.0, 0.0, 0.0});
        MathLib::setMatrix(K, {1.0, 0.0, 0.0, 1.0});
        MathLib::setVector(b, {x0 - 0.5 * x1 + 1.0, x0});
    }

    bool isLinear() const override { return false; }
};
}  // namespace

/*! The delta-x convergence must be checked on the Picard step, not on the
 * mixed one: here the mixed step of the second iteration vanishes, and a check
 * after the mixing reports convergence at x = (1, 0).
 */
// Disabled for PETSc together with the other ODE system tests, see issue
// #1989.
#ifndef USE_PETSC
TEST(NumLibNonlinearSolver, PicardAndersonVanishingMixedStepIsNotConvergence)
#else
TEST(NumLibNonlinearSolver,
     DISABLED_PicardAndersonVanishingMixedStepIsNotConvergence)
#endif
{
    int const process_id = 0;
    int const maxiter = 50;
    int const anderson_depth = 2;

    ODEMixingReproducesIterate ode;
    NumLib::BackwardEuler time_disc;
    NumLib::TimeDiscretizedODESystem<ODEMixingReproducesIterate::ODETag,
                                     NumLib::NonlinearSolverTag::Picard>
        ode_sys(process_id, ode, time_disc);

    auto const linear_solver = createLinearSolver();
    NumLib::ConvergenceCriterionDeltaX convergence_criterion(
        1e-9, std::nullopt, MathLib::VecNormType::NORM2);
    NumLib::NonlinearSolver<NumLib::NonlinearSolverTag::Picard>
        nonlinear_solver(*linear_solver, maxiter, anderson_depth, 1.0);
    nonlinear_solver.setEquationSystem(ode_sys, convergence_criterion);

    GlobalVector x0(ode.getMatrixSpecifications(process_id).nrows);
    MathLib::setVector(x0, {0.0, 0.0});
    MathLib::LinAlg::finalizeAssembly(x0);

    auto& x = NumLib::GlobalVectorProvider::provider.getVector(x0);
    auto& x_prev = NumLib::GlobalVectorProvider::provider.getVector(x0);
    std::vector<GlobalVector*> xs{&x};
    std::vector<GlobalVector*> xs_prev{&x_prev};

    time_disc.setInitialState(0.0);
    time_disc.nextTimestep(1.0, 1.0);

    auto const status =
        nonlinear_solver.solve(xs, xs_prev, nullptr, process_id);

    EXPECT_TRUE(status.error_norms_met);
    EXPECT_GT(status.number_iterations, 2);
    MathLib::LinAlg::setLocalAccessibleVector(x);
    EXPECT_NEAR(x.get(0), 2.0, 1e-8);
    EXPECT_NEAR(x.get(1), 2.0, 1e-8);

    NumLib::GlobalVectorProvider::provider.releaseVector(x);
    NumLib::GlobalVectorProvider::provider.releaseVector(x_prev);
}
