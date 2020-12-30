/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <functional>
#include <memory>


#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"
#include "ProcessLib/Output/Output.h"

#include "Process.h"

namespace NumLib
{
class ConvergenceCriterion;
}

namespace ChemistryLib
{
class ChemicalSolverInterface;
}

namespace ProcessLib
{
struct ProcessData;

/// Time loop capable of time-integrating several processes at once.
class TimeLoop
{
public:
    TimeLoop(std::unique_ptr<Output>&& output,
             std::vector<std::unique_ptr<ProcessData>>&& per_process_data,
             const int global_coupling_max_iterations,
             std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>&&
                 global_coupling_conv_crit,
             const double start_time, const double end_time);

    void initialize();
    bool loop();

    ~TimeLoop();

private:
    /**
     * This function fills the vector of solutions of coupled processes of
     * processes, _solutions_of_coupled_processes, and initializes the vector
     * of solutions of the previous coupling iteration,
     * _solutions_of_last_cpl_iteration.
     */
    void setCoupledSolutions();

    /**
     * \brief Member to solver non coupled systems of equations, which can be
     *        a single system of equations, or several systems of equations
     *        without any dependency among the different systems.
     *
     * @param t           Current time
     * @param dt          Time step size
     * @param timestep_id Index of the time step
     * @return            true:  if all nonlinear solvers convergence.
     *                    false: if any of nonlinear solvers divergences.
     */
    NumLib::NonlinearSolverStatus solveUncoupledEquationSystems(
        const double t, const double dt, const std::size_t timestep_id);

    /**
     * \brief Member to solver coupled systems of equations by the staggered
     *        scheme.
     *
     * @param t           Current time
     * @param dt          Time step size
     * @param timestep_id Index of the time step
     * @return            true:   if all nonlinear solvers convergence.
     *                    false:  if any of nonlinear solvers divergences.
     */
    NumLib::NonlinearSolverStatus solveCoupledEquationSystemsByStaggeredScheme(
        const double t, const double dt, const std::size_t timestep_id);

    /**
     *  Find the minimum time step size among the predicted step sizes of
     *  processes and step it as common time step size.
     *
     *  @param prev_dt        Previous time step size.
     *  @param t              Current time.
     *  @param accepted_steps Accepted time steps that are counted in this
     *                        function.
     *  @param rejected_steps Rejected time steps that are counted in this
     *                        function.
     */
    double computeTimeStepping(const double prev_dt, double& t,
                               std::size_t& accepted_steps,
                               std::size_t& rejected_steps);

    template <typename OutputClass, typename OutputClassMember>
    void outputSolutions(bool const output_initial_condition, unsigned timestep,
                         const double t, OutputClass& output_object,
                         OutputClassMember output_class_member) const;

private:
    std::vector<GlobalVector*> _process_solutions;
    std::vector<GlobalVector*> _process_solutions_prev;
    std::unique_ptr<Output> _output;
    std::vector<std::unique_ptr<ProcessData>> _per_process_data;

    bool _last_step_rejected = false;
    int _repeating_times_of_rejected_step = 0;
    const double _start_time;
    const double _end_time;

    /// Maximum iterations of the global coupling.
    const int _global_coupling_max_iterations;
    /// Convergence criteria of processes for the global coupling iterations.
    std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>
        _global_coupling_conv_crit;

    /// Solutions of the previous coupling iteration for the convergence
    /// criteria of the coupling iteration.
    std::vector<GlobalVector*> _solutions_of_last_cpl_iteration;
};
}  // namespace ProcessLib
