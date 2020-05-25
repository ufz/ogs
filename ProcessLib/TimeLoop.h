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
             std::unique_ptr<ChemistryLib::ChemicalSolverInterface>&&
                 chemical_system,
             const double start_time, const double end_time);

    void initialize();
    bool loop();

    ~TimeLoop();

private:
    /**
     * This function fills the vector of solutions of coupled processes of
     * processes, solutions_of_coupled_processes_, and initializes the vector
     * of solutions of the previous coupling iteration,
     * solutions_of_last_cpl_iteration_.
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
    std::vector<GlobalVector*> process_solutions_;
    std::vector<GlobalVector*> process_solutions_prev_;
    std::unique_ptr<Output> output_;
    std::vector<std::unique_ptr<ProcessData>> per_process_data_;

    bool last_step_rejected_ = false;
    int repeating_times_of_rejected_step_ = 0;
    const double start_time_;
    const double end_time_;

    /// Maximum iterations of the global coupling.
    const int global_coupling_max_iterations_;
    /// Convergence criteria of processes for the global coupling iterations.
    std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>
        global_coupling_conv_crit_;

    std::unique_ptr<ChemistryLib::ChemicalSolverInterface>
        chemical_solver_interface_;

    /// Solutions of the previous coupling iteration for the convergence
    /// criteria of the coupling iteration.
    std::vector<GlobalVector*> solutions_of_last_cpl_iteration_;
};
}  // namespace ProcessLib
