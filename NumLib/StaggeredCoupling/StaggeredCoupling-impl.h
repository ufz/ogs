/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 20, 2023, 5:09 PM
 */
#pragma once

#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "NumLib/ODESolver/NonlinearSolver.h"
#include "StaggeredCoupling.h"

namespace NumLib
{
template <typename ProcessData, typename Output>
NumLib::NonlinearSolverStatus StaggeredCoupling::execute(
    const double t, const double dt, const std::size_t timestep_id,
    std::vector<GlobalVector*>& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<Output> const& outputs,
    ProcessSolver<ProcessData, Output> const& solve_one_time_step_one_process)
{
    auto const [nonlinear_solver_status, coupling_iteration_converged] =
        executeConcrete(coupling_nodes_, global_coupling_max_iterations_, t, dt,
                        timestep_id, process_solutions, process_solutions_prev,
                        per_process_data, outputs,
                        solve_one_time_step_one_process);

    if (!coupling_iteration_converged)
    {
        WARN(
            "The coupling iterations reaches its maximum number in time step "
            "#{:d} at t = {:g} s",
            timestep_id, t);
    }

    return nonlinear_solver_status;
}

template <typename ProcessData, typename Output>
NumLib::NonlinearSolverStatus StaggeredCoupling::executeSingleIteration(
    int const global_coupling_iteration,
    CouplingNode const& regular_coupling_node, const double t, const double dt,
    const std::size_t timestep_id,
    std::vector<GlobalVector*>& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<Output> const& outputs,
    ProcessSolver<ProcessData, Output> const& solve_one_time_step_one_process)
{
    BaseLib::RunTime time_timestep_process;
    time_timestep_process.start();

    auto const process_id = regular_coupling_node.process_id;
    auto const& process_name = regular_coupling_node.process_name;
    INFO("Solve process #{:d} (named as {:s})", process_id, process_name);

    auto& process_data = *(per_process_data[process_id]);

    process_data.nonlinear_solver_status = solve_one_time_step_one_process(
        process_solutions, process_solutions_prev, timestep_id, t, dt,
        process_data, outputs);

    INFO(
        "[time] Solving process #{:d} (named as {:s}) took {:g} s in "
        "time step #{} coupling iteration #{}.",
        process_id, process_name, time_timestep_process.elapsed(), timestep_id,
        global_coupling_iteration);

    return process_data.nonlinear_solver_status;
}

template <typename ProcessData, typename Output>
std::tuple<NumLib::NonlinearSolverStatus, bool>
StaggeredCoupling::executeConcrete(
    std::vector<CouplingNodeVariant>& coupling_nodes, const int max_iterations,
    const double t, const double dt, const std::size_t timestep_id,
    std::vector<GlobalVector*>& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<Output> const& outputs,
    ProcessSolver<ProcessData, Output> const& solve_one_time_step_one_process)
{
    setFirstIterationIndicator(coupling_nodes);

    NumLib::NonlinearSolverStatus nonlinear_solver_status{true, -1};

    bool coupling_iteration_converged = true;
    for (int global_coupling_iteration = 0;
         global_coupling_iteration < max_iterations;
         global_coupling_iteration++,
             resetCouplingConvergenceCriteria(coupling_nodes))
    {
        coupling_iteration_converged = true;
        for (auto& coupling_node : coupling_nodes)
        {
            // For the dummy root node, perform sub-coupling computation.
            if (std::holds_alternative<RootCouplingNode>(coupling_node))
            {
                auto const [local_nonlinear_solver_status,
                            local_coupling_iteration_converged] =
                    executeSubCoupling(
                        coupling_node, t, dt, timestep_id, process_solutions,
                        process_solutions_prev, per_process_data, outputs,
                        solve_one_time_step_one_process);

                if (!local_nonlinear_solver_status.error_norms_met)
                {
                    coupling_iteration_converged = false;
                    return {local_nonlinear_solver_status,
                            coupling_iteration_converged};
                }

                coupling_iteration_converged =
                    coupling_iteration_converged &&
                    local_coupling_iteration_converged;
                continue;
            }

            CouplingNode const& regular_coupling_node =
                std::get<CouplingNode>(coupling_node);

            nonlinear_solver_status = executeSingleIteration(
                global_coupling_iteration, regular_coupling_node, t, dt,
                timestep_id, process_solutions, process_solutions_prev,
                per_process_data, outputs, solve_one_time_step_one_process);

            if (!nonlinear_solver_status.error_norms_met)
            {
                WARN(
                    "The nonlinear solver failed in time step #{:d} at t = "
                    "{:g} s for process {:s}.",
                    timestep_id, t, regular_coupling_node.process_name);
                coupling_iteration_converged = false;
                return {nonlinear_solver_status, coupling_iteration_converged};
            }

            auto const& x =
                *process_solutions[regular_coupling_node.process_id];

            // It is unnecessary to check convergence for a process at the first
            // iteration
            if (global_coupling_iteration > 0)
            {
                coupling_iteration_converged = checkCouplingConvergence(
                    coupling_iteration_converged, regular_coupling_node, x);
            }

            updatePreviousSolution(regular_coupling_node.process_id, x);

        }  // end of for (auto& process_data : _per_process_data)

        // At least to run two coupling iterations, meaning that the coupling
        // has at least two coupling nodes.
        if (coupling_iteration_converged && global_coupling_iteration > 0)
        {
            break;
        }
    }

    return {nonlinear_solver_status, coupling_iteration_converged};
}

template <typename ProcessData, typename Output>
std::tuple<NumLib::NonlinearSolverStatus, bool>
StaggeredCoupling::executeSubCoupling(
    CouplingNodeVariant& coupling_node, const double t, const double dt,
    const std::size_t timestep_id,
    std::vector<GlobalVector*>& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<Output> const& outputs,
    ProcessSolver<ProcessData, Output> const& solve_one_time_step_one_process)
{
    INFO("--- Execute sub-coupling:");
    RootCouplingNode& root_coupling_node =
        std::get<RootCouplingNode>(coupling_node);
    const int local_max_iterations =
        std::get<CouplingNode>(root_coupling_node.sub_coupling_nodes.front())
            .max_iterations;

    auto const [sub_nonlinear_solver_status, sub_coupling_iteration_converged] =
        executeConcrete<ProcessData, Output>(
            root_coupling_node.sub_coupling_nodes, local_max_iterations, t, dt,
            timestep_id, process_solutions, process_solutions_prev,
            per_process_data, outputs, solve_one_time_step_one_process);

    INFO("--- End sub-coupling.");
    return {sub_nonlinear_solver_status, sub_coupling_iteration_converged};
}

}  // namespace NumLib
