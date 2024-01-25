/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <functional>
#include <memory>

#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"
#include "Process.h"
#include "ProcessLib/Output/Output.h"

namespace NumLib
{
class ConvergenceCriterion;
class StaggeredCoupling;
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
    TimeLoop(std::vector<Output>&& outputs,
             std::vector<std::unique_ptr<ProcessData>>&& per_process_data,
             std::unique_ptr<NumLib::StaggeredCoupling>&& staggered_coupling,
             const double start_time, const double end_time);

    void initialize();
    void outputLastTimeStep() const;

    ~TimeLoop();

    bool executeTimeStep();

    /// Computes and sets the next timestep.
    ///
    /// \attention The timestepper might reject the current timestep and repeat
    /// it (with a reduced timestep size).
    ///
    /// \returns true if the simulation (time) has not finished, yet, false
    /// otherwise.
    bool calculateNextTimeStep();

    double endTime() const { return _end_time; }
    double currentTime() const { return _current_time; }
    bool successful_time_step = false;

private:
    bool preTsNonlinearSolvePostTs(double const t, double const dt,
                                   std::size_t const timesteps);

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
     *  @param time_step_constraints Functions that are evaluate to
     *  influence the time step size (for instance a fixed output time)
     *  @return the time step size and the information if the last time step was
     *  rejected
     */
    std::pair<double, bool> computeTimeStepping(
        const double prev_dt, double& t, std::size_t& accepted_steps,
        std::size_t& rejected_steps,
        std::vector<std::function<double(double, double)>> const&
            time_step_constraints);

    template <typename OutputClassMember>
    void outputSolutions(unsigned timestep,
                         const double t,
                         OutputClassMember output_class_member) const;

private:
    std::vector<std::function<double(double, double)>>
    generateOutputTimeStepConstraints(std::vector<double>&& fixed_times) const;
    void preOutputInitialConditions(const double t) const;
    std::vector<GlobalVector*> _process_solutions;
    std::vector<GlobalVector*> _process_solutions_prev;
    std::vector<Output> _outputs;
    std::vector<std::unique_ptr<ProcessData>> _per_process_data;

    const double _start_time;
    const double _end_time;
    double _current_time = _start_time;
    std::size_t _accepted_steps = 0;
    std::size_t _rejected_steps = 0;
    double _dt = 0;
    int _repeating_times_of_rejected_step = 0;
    bool _last_step_rejected = false;

    std::unique_ptr<NumLib::StaggeredCoupling> _staggered_coupling;
};
}  // namespace ProcessLib
