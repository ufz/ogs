// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/ODESolver/ConvergenceCriterion.h"
#include "NumLib/ODESolver/NonlinearSolverStatus.h"
#include "NumLib/StaggeredCoupling/StaggeredCoupling.h"

/// Building blocks for driving NumLib::StaggeredCoupling without any real
/// process: convergence criteria with a prescribed outcome, process data,
/// process solvers, and the coupled system they belong to.
namespace StaggeredCouplingTestDoubles
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

/// Two coupling nodes, one per process, each allowed one iteration.
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

inline std::vector<std::unique_ptr<TestProcessData>> createProcessData()
{
    std::vector<std::unique_ptr<TestProcessData>> per_process_data;
    per_process_data.emplace_back(std::make_unique<TestProcessData>());
    per_process_data.emplace_back(std::make_unique<TestProcessData>());
    return per_process_data;
}

/// Solves nothing; it only reports the given status and changes the solution
/// so that the coupling convergence check operates on a non-zero increment.
inline auto createProcessSolver(bool const nonlinear_solver_converged)
{
    return [nonlinear_solver_converged, value = 0.0](
               std::vector<GlobalVector*> const& xs,
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

/// The two coupled processes a NumLib::StaggeredCoupling is executed over:
/// their solution vectors, their process data and their outputs.
///
/// Everything NumLib::StaggeredCoupling::execute() needs from its caller, and
/// nothing else -- in particular no time stepping. Whoever the caller is, be
/// it ProcessLib::TimeLoop or this test, the coupling sees the same
/// arguments.
struct CoupledSystem
{
    CoupledSystem()
    {
        for (auto* x : {&x0, &x1, &x0_prev, &x1_prev})
        {
            MathLib::LinAlg::set(*x, 0.0);
        }
    }

    /// Runs the coupling scheme once, over a solve which reports the given
    /// convergence status.
    NumLib::NonlinearSolverStatus solveOneStep(
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
}  // namespace StaggeredCouplingTestDoubles
