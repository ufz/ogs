/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TimeLoop.h"

#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "ChemistryLib/CreatePhreeqcIO.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "ProcessLib/CreateProcessData.h"
#include "ProcessLib/Output/CreateOutput.h"

#include "CoupledSolutionsForStaggeredScheme.h"
#include "ProcessData.h"

namespace
{
//! Sets the EquationSystem for the given nonlinear solver,
//! which is Picard or Newton depending on the NLTag.
template <NumLib::NonlinearSolverTag NLTag>
void setEquationSystem(NumLib::NonlinearSolverBase& nonlinear_solver,
                              NumLib::EquationSystem& eq_sys,
                              NumLib::ConvergenceCriterion& conv_crit)
{
    using Solver = NumLib::NonlinearSolver<NLTag>;
    using EqSys = NumLib::NonlinearSystem<NLTag>;

    assert(dynamic_cast<Solver*>(&nonlinear_solver) != nullptr);
    assert(dynamic_cast<EqSys*>(&eq_sys) != nullptr);

    auto& nl_solver_ = static_cast<Solver&>(nonlinear_solver);
    auto& eq_sys_ = static_cast<EqSys&>(eq_sys);

    nl_solver_.setEquationSystem(eq_sys_, conv_crit);
}

//! Sets the EquationSystem for the given nonlinear solver,
//! transparently both for Picard and Newton solvers.
void setEquationSystem(NumLib::NonlinearSolverBase& nonlinear_solver,
                              NumLib::EquationSystem& eq_sys,
                              NumLib::ConvergenceCriterion& conv_crit,
                              NumLib::NonlinearSolverTag nl_tag)
{
    using Tag = NumLib::NonlinearSolverTag;
    switch (nl_tag)
    {
        case Tag::Picard:
            setEquationSystem<Tag::Picard>(nonlinear_solver, eq_sys, conv_crit);
            break;
        case Tag::Newton:
            setEquationSystem<Tag::Newton>(nonlinear_solver, eq_sys, conv_crit);
            break;
    }
}

}  // namespace


namespace ProcessLib
{
template <NumLib::ODESystemTag ODETag>
void setTimeDiscretizedODESystem(
    ProcessData& process_data,
    NumLib::ODESystem<ODETag, NumLib::NonlinearSolverTag::Picard>& ode_sys)
{
    using Tag = NumLib::NonlinearSolverTag;
    // A concrete Picard solver
    using NonlinearSolverPicard = NumLib::NonlinearSolver<Tag::Picard>;
    // A concrete Newton solver
    using NonlinearSolverNewton = NumLib::NonlinearSolver<Tag::Newton>;

    if (dynamic_cast<NonlinearSolverPicard*>(&process_data.nonlinear_solver))
    {
        // The Picard solver can also work with a Newton-ready ODE,
        // because the Newton ODESystem derives from the Picard ODESystem.
        // So no further checks are needed here.

        process_data.tdisc_ode_sys = std::make_unique<
            NumLib::TimeDiscretizedODESystem<ODETag, Tag::Picard>>(
            process_data.process_id, ode_sys, *process_data.time_disc);
    }
    else if (dynamic_cast<NonlinearSolverNewton*>(
                 &process_data.nonlinear_solver))
    {
        // The Newton-Raphson method needs a Newton-ready ODE.

        using ODENewton = NumLib::ODESystem<ODETag, Tag::Newton>;
        if (auto* ode_newton = dynamic_cast<ODENewton*>(&ode_sys))
        {
            process_data.tdisc_ode_sys = std::make_unique<
                NumLib::TimeDiscretizedODESystem<ODETag, Tag::Newton>>(
                process_data.process_id, *ode_newton, *process_data.time_disc);
        }
        else
        {
            OGS_FATAL(
                "You are trying to solve a non-Newton-ready ODE with the"
                " Newton-Raphson method. Aborting");
        }
    }
    else
    {
        OGS_FATAL("Encountered unknown nonlinear solver type. Aborting");
    }

    process_data.mat_strg = dynamic_cast<NumLib::InternalMatrixStorage*>(
        process_data.tdisc_ode_sys.get());
}

void setTimeDiscretizedODESystem(ProcessData& process_data)
{
    setTimeDiscretizedODESystem(process_data, process_data.process);
}

std::vector<GlobalVector*> setInitialConditions(
    double const t0,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data)
{
    std::vector<GlobalVector*> process_solutions;

    int process_id = 0;
    for (auto& process_data : per_process_data)
    {
        auto& pcs = process_data->process;
        auto& time_disc = *process_data->time_disc;

        auto& ode_sys = *process_data->tdisc_ode_sys;
        auto const nl_tag = process_data->nonlinear_solver_tag;

        // append a solution vector of suitable size
        process_solutions.emplace_back(
            &NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications(process_id)));

        auto& x0 = *process_solutions[process_id];
        pcs.setInitialConditions(process_id, t0, x0);
        MathLib::LinAlg::finalizeAssembly(x0);

        time_disc.setInitialState(t0, x0);  // push IC

        if (time_disc.needsPreload())
        {
            auto& nonlinear_solver = process_data->nonlinear_solver;
            auto& mat_strg = *process_data->mat_strg;
            auto& conv_crit = *process_data->conv_crit;

            setEquationSystem(nonlinear_solver, ode_sys, conv_crit, nl_tag);
            nonlinear_solver.assemble(x0);
            time_disc.pushState(
                t0, x0,
                mat_strg);  // TODO: that might do duplicate work
        }

        ++process_id;
    }

    return process_solutions;
}

NumLib::NonlinearSolverStatus solveOneTimeStepOneProcess(
    int const process_id, GlobalVector& x, std::size_t const timestep,
    double const t, double const delta_t, ProcessData const& process_data,
    Output& output_control)
{
    auto& process = process_data.process;
    auto& time_disc = *process_data.time_disc;
    auto& conv_crit = *process_data.conv_crit;
    auto& ode_sys = *process_data.tdisc_ode_sys;
    auto& nonlinear_solver = process_data.nonlinear_solver;
    auto const nl_tag = process_data.nonlinear_solver_tag;

    setEquationSystem(nonlinear_solver, ode_sys, conv_crit, nl_tag);

    // Note: Order matters!
    // First advance to the next timestep, then set known solutions at that
    // time, afterwards pass the right solution vector and time to the
    // preTimestep() hook.

    time_disc.nextTimestep(t, delta_t);

    auto const post_iteration_callback = [&](int iteration,
                                             GlobalVector const& x) {
        output_control.doOutputNonlinearIteration(process, process_id, timestep,
                                                  t, x, iteration);
    };

    auto const nonlinear_solver_status =
        nonlinear_solver.solve(x, post_iteration_callback);

    if (nonlinear_solver_status.error_norms_met)
    {
        process.postNonLinearSolver(x, t, process_id);
    }

    return nonlinear_solver_status;
}

TimeLoop::TimeLoop(std::unique_ptr<Output>&& output,
                   std::vector<std::unique_ptr<ProcessData>>&& per_process_data,
                   const int global_coupling_max_iterations,
                   std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>&&
                       global_coupling_conv_crit,
                   std::unique_ptr<ChemistryLib::PhreeqcIO>&& chemical_system,
                   const double start_time, const double end_time)
    : _output(std::move(output)),
      _per_process_data(std::move(per_process_data)),
      _start_time(start_time),
      _end_time(end_time),
      _global_coupling_max_iterations(global_coupling_max_iterations),
      _global_coupling_conv_crit(std::move(global_coupling_conv_crit)),
      _chemical_system(std::move(chemical_system))
{
}

bool TimeLoop::setCoupledSolutions()
{
    // All _per_process_data share one process
    const bool use_monolithic_scheme =
        _per_process_data[0]->process.isMonolithicSchemeUsed();
    if (use_monolithic_scheme)
    {
        return false;
    }

    _solutions_of_coupled_processes.reserve(_per_process_data.size());
    for (unsigned process_id = 0; process_id < _per_process_data.size();
         process_id++)
    {
        auto const& x = *_process_solutions[process_id];
        _solutions_of_coupled_processes.emplace_back(x);

        // Create a vector to store the solution of the last coupling iteration
        auto& x0 = NumLib::GlobalVectorProvider::provider.getVector(x);
        MathLib::LinAlg::copy(x, x0);

        // append a solution vector of suitable size
        _solutions_of_last_cpl_iteration.emplace_back(&x0);
    }

    return true;  // use staggered scheme.
}

double TimeLoop::computeTimeStepping(const double prev_dt, double& t,
                                     std::size_t& accepted_steps,
                                     std::size_t& rejected_steps)
{
    bool all_process_steps_accepted = true;
    // Get minimum time step size among step sizes of all processes.
    double dt = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < _per_process_data.size(); i++)
    {
        auto& ppd = *_per_process_data[i];
        const auto& timestepper = ppd.timestepper;

        auto& time_disc = ppd.time_disc;
        auto const& x = *_process_solutions[i];

        auto const& conv_crit = ppd.conv_crit;
        const MathLib::VecNormType norm_type =
            (conv_crit) ? conv_crit->getVectorNormType()
                        : MathLib::VecNormType::NORM2;

        const double solution_error =
            (timestepper->isSolutionErrorComputationNeeded())
                ? ((t == timestepper->begin())
                       ? 0.  // Always accepts the zeroth step
                       : time_disc->getRelativeChangeFromPreviousTimestep(
                             x, norm_type))
                : 0.;

        if (!ppd.nonlinear_solver_status.error_norms_met)
        {
            timestepper->setAcceptedOrNot(false);
        }
        else
        {
            timestepper->setAcceptedOrNot(true);
        }

        if (!timestepper->next(solution_error,
                               ppd.nonlinear_solver_status.number_iterations) &&
            // In case of FixedTimeStepping, which makes timestepper->next(...)
            // return false when the ending time is reached.
            t + std::numeric_limits<double>::epsilon() < timestepper->end())
        {
            // Not all processes have accepted steps.
            all_process_steps_accepted = false;
        }

        if (!ppd.nonlinear_solver_status.error_norms_met)
        {
            WARN("Time step will be rejected due to nonlinear solver diverged");
            all_process_steps_accepted = false;
        }

        if (timestepper->getTimeStep().dt() >
                std::numeric_limits<double>::min() ||
            std::abs(t - timestepper->end()) <
                std::numeric_limits<double>::epsilon())
        {
            if (timestepper->getTimeStep().dt() < dt)
            {
                dt = timestepper->getTimeStep().dt();
            }
        }
        else
        {
            // dt being close to 0 only happens when
            // t_n + dt > t_s, and dt is forced to be zero. Where t_n the time
            // of previous time step, and t_s is the specified time taken from
            // input or the end time. Under this condition, the time stepping
            // is skipped.
            ppd.skip_time_stepping = true;
        }
    }

    if (all_process_steps_accepted)
    {
        _repeating_times_of_rejected_step = 0;
    }
    else
    {
        _repeating_times_of_rejected_step++;
    }

    bool is_initial_step = false;
    // Reset the time step with the minimum step size, dt
    // Update the solution of the previous time step in time_disc.
    for (std::size_t i = 0; i < _per_process_data.size(); i++)
    {
        const auto& ppd = *_per_process_data[i];
        auto& timestepper = ppd.timestepper;
        timestepper->resetCurrentTimeStep(dt);

        if (ppd.skip_time_stepping)
        {
            continue;
        }

        if (t == timestepper->begin())
        {
            is_initial_step = true;
            continue;
        }

        auto& time_disc = ppd.time_disc;
        auto& x = *_process_solutions[i];
        if (all_process_steps_accepted)
        {
            time_disc->pushState(t, x, *ppd.mat_strg);
        }
        else
        {
            if (t < _end_time || std::abs(t - _end_time) <
                                     std::numeric_limits<double>::epsilon())
            {
                WARN(
                    "Time step %d was rejected %d times "
                    "and it will be repeated with a reduced step size.",
                    accepted_steps + 1, _repeating_times_of_rejected_step);
                time_disc->popState(x);
            }
        }
    }

    if (!is_initial_step)
    {
        if (all_process_steps_accepted)
        {
            accepted_steps++;
            _last_step_rejected = false;
        }
        else
        {
            if (t < _end_time || std::abs(t - _end_time) <
                                     std::numeric_limits<double>::epsilon())
            {
                t -= prev_dt;
                rejected_steps++;
                _last_step_rejected = true;
            }
        }
    }

    // Adjust step size if t < _end_time, while t+dt exceeds the end time
    if (t < _end_time && t + dt > _end_time)
    {
        dt = _end_time - t;
    }

    return dt;
}

/*
 * TODO:
 * Now we have a structure inside the time loop which is very similar to the
 * nonlinear solver. And admittedly, the control flow inside the nonlinear
 * solver is rather complicated. Maybe in the future con can introduce an
 * abstraction that can do both the convergence checks of the coupling loop and
 * of the nonlinear solver.
 *
 */
bool TimeLoop::loop()
{
    // initialize output, convergence criterion, etc.
    {
        int process_id = 0;
        for (auto& process_data : _per_process_data)
        {
            auto& pcs = process_data->process;
            _output->addProcess(pcs, process_id);

            process_data->process_id = process_id;
            setTimeDiscretizedODESystem(*process_data);

            if (auto* conv_crit =
                    dynamic_cast<NumLib::ConvergenceCriterionPerComponent*>(
                        process_data->conv_crit.get()))
            {
                conv_crit->setDOFTable(pcs.getDOFTable(process_id),
                                       pcs.getMesh());
            }

            // Add the fixed times of output to time stepper in order that
            // the time stepping is performed and the results are output at
            // these times. Note: only the adaptive time steppers can have the
            // the fixed times.
            auto& timestepper = process_data->timestepper;
            timestepper->addFixedOutputTimes(_output->getFixedOutputTimes());

            ++process_id;
        }
    }

    // init solution storage
    _process_solutions = setInitialConditions(_start_time, _per_process_data);

    const bool is_staggered_coupling = setCoupledSolutions();

    // Output initial conditions
    {
        const bool output_initial_condition = true;
        outputSolutions(output_initial_condition, is_staggered_coupling, 0,
                        _start_time, *_output, &Output::doOutput);
    }

    double t = _start_time;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    NumLib::NonlinearSolverStatus nonlinear_solver_status{true, 0};

    double dt = computeTimeStepping(0.0, t, accepted_steps, rejected_steps);

    while (t < _end_time)
    {
        BaseLib::RunTime time_timestep;
        time_timestep.start();

        t += dt;
        const double prev_dt = dt;

        const std::size_t timesteps = accepted_steps + 1;
        // TODO(wenqing): , input option for time unit.
        INFO("=== Time stepping at step #%u and time %g with step size %g",
             timesteps, t, dt);

        // Check element deactivation:
        int process_id = 0;
        for (auto& process_data : _per_process_data)
        {
            process_data->process.updateDeactivatedSubdomains(t, process_id);
            ++process_id;
        }

        if (is_staggered_coupling)
        {
            nonlinear_solver_status =
                solveCoupledEquationSystemsByStaggeredScheme(t, dt, timesteps);
        }
        else
        {
            nonlinear_solver_status =
                solveUncoupledEquationSystems(t, dt, timesteps);
        }

        INFO("[time] Time step #%u took %g s.", timesteps,
             time_timestep.elapsed());

        dt = computeTimeStepping(prev_dt, t, accepted_steps, rejected_steps);

        if (!_last_step_rejected)
        {
            const bool output_initial_condition = false;
            outputSolutions(output_initial_condition, is_staggered_coupling,
                            timesteps, t, *_output, &Output::doOutput);
        }

        if (t + dt > _end_time ||
            t + std::numeric_limits<double>::epsilon() > _end_time)
        {
            break;
        }

        if (dt < std::numeric_limits<double>::epsilon())
        {
            WARN(
                "Time step size of %g is too small.\n"
                "Time stepping stops at step %u and at time of %g.",
                dt, timesteps, t);
            break;
        }
    }

    INFO(
        "The whole computation of the time stepping took %u steps, in which\n"
        "\t the accepted steps are %u, and the rejected steps are %u.\n",
        accepted_steps + rejected_steps, accepted_steps, rejected_steps);

    // output last time step
    if (nonlinear_solver_status.error_norms_met)
    {
        const bool output_initial_condition = false;
        outputSolutions(output_initial_condition, is_staggered_coupling,
                        accepted_steps + rejected_steps, t, *_output,
                        &Output::doOutputLastTimestep);
    }

    return nonlinear_solver_status.error_norms_met;
}

static std::string const nonlinear_fixed_dt_fails_info =
    "Nonlinear solver fails. Because the time stepper FixedTimeStepping is "
    "used, the program has to be terminated.";

NumLib::NonlinearSolverStatus TimeLoop::solveUncoupledEquationSystems(
    const double t, const double dt, const std::size_t timestep_id)
{
    NumLib::NonlinearSolverStatus nonlinear_solver_status;
    // TODO(wenqing): use process name
    unsigned process_id = 0;
    for (auto& process_data : _per_process_data)
    {
        if (process_data->skip_time_stepping)
        {
            INFO("Process %u is skipped in the time stepping.", process_id);
            ++process_id;
            continue;
        }

        BaseLib::RunTime time_timestep_process;
        time_timestep_process.start();

        auto& x = *_process_solutions[process_id];
        auto& pcs = process_data->process;
        pcs.preTimestep(x, t, dt, process_id);

        nonlinear_solver_status = solveOneTimeStepOneProcess(
            process_id, x, timestep_id, t, dt, *process_data, *_output);
        process_data->nonlinear_solver_status = nonlinear_solver_status;
        pcs.postTimestep(x, t, dt, process_id);
        pcs.computeSecondaryVariable(t, x, process_id);

        INFO("[time] Solving process #%u took %g s in time step #%u ",
             process_id, time_timestep_process.elapsed(), timestep_id);

        if (!nonlinear_solver_status.error_norms_met)
        {
            ERR("The nonlinear solver failed in time step #%u at t = %g s for "
                "process #%u.",
                timestep_id, t, process_id);

            if (!process_data->timestepper->isSolutionErrorComputationNeeded())
            {
                // save unsuccessful solution
                _output->doOutputAlways(pcs, process_id, timestep_id, t, x);
                OGS_FATAL(nonlinear_fixed_dt_fails_info.data());
            }

            return nonlinear_solver_status;
        }

        ++process_id;
    }  // end of for (auto& process_data : _per_process_data)

    return nonlinear_solver_status;
}

NumLib::NonlinearSolverStatus
TimeLoop::solveCoupledEquationSystemsByStaggeredScheme(
    const double t, const double dt, const std::size_t timestep_id)
{
    // Coupling iteration
    if (_global_coupling_max_iterations != 0)
    {
        // Set the flag of the first iteration be true.
        for (auto& conv_crit : _global_coupling_conv_crit)
        {
            conv_crit->preFirstIteration();
        }
    }
    auto resetCouplingConvergenceCriteria = [&]() {
        for (auto& conv_crit : _global_coupling_conv_crit)
        {
            conv_crit->reset();
        }
    };

    // Update solutions of previous time step at once
    {
        int process_id = 0;
        for (auto& process_data : _per_process_data)
        {
            auto& x = *_process_solutions[process_id];
            process_data->process.preTimestep(x, t, dt, process_id);
            ++process_id;
        }
    }

    NumLib::NonlinearSolverStatus nonlinear_solver_status{true, 0};
    bool coupling_iteration_converged = true;
    for (int global_coupling_iteration = 0;
         global_coupling_iteration < _global_coupling_max_iterations;
         global_coupling_iteration++, resetCouplingConvergenceCriteria())
    {
        // TODO(wenqing): use process name
        coupling_iteration_converged = true;
        int process_id = 0;
        int const last_process_id = _per_process_data.size() - 1;
        for (auto& process_data : _per_process_data)
        {
            if (process_data->skip_time_stepping)
            {
                INFO("Process %u is skipped in the time stepping.", process_id);
                ++process_id;
                continue;
            }

            BaseLib::RunTime time_timestep_process;
            time_timestep_process.start();

            auto& x = *_process_solutions[process_id];

            CoupledSolutionsForStaggeredScheme coupled_solutions(
                _solutions_of_coupled_processes, dt, process_id);

            process_data->process.setCoupledSolutionsForStaggeredScheme(
                &coupled_solutions);

            nonlinear_solver_status = solveOneTimeStepOneProcess(
                process_id, x, timestep_id, t, dt, *process_data, *_output);
            process_data->nonlinear_solver_status = nonlinear_solver_status;

            INFO(
                "[time] Solving process #%u took %g s in time step #%u "
                " coupling iteration #%u",
                process_id, time_timestep_process.elapsed(), timestep_id,
                global_coupling_iteration);

            if (!nonlinear_solver_status.error_norms_met)
            {
                ERR("The nonlinear solver failed in time step #%u at t = %g s "
                    "for process #%u.",
                    timestep_id, t, process_id);

                if (!process_data->timestepper
                         ->isSolutionErrorComputationNeeded())
                {
                    // save unsuccessful solution
                    _output->doOutputAlways(process_data->process, process_id,
                                            timestep_id, t, x);
                    OGS_FATAL(nonlinear_fixed_dt_fails_info.data());
                }
                break;
            }

            // Check the convergence of the coupling iteration
            auto& x_old = *_solutions_of_last_cpl_iteration[process_id];
            if (global_coupling_iteration > 0)
            {
                MathLib::LinAlg::axpy(x_old, -1.0, x);  // save dx to x_old
                if (process_id == last_process_id)
                {
                    INFO(
                        "------- Checking convergence criterion for coupled "
                        "solution  -------");
                    _global_coupling_conv_crit[process_id]->checkDeltaX(x_old,
                                                                        x);
                    coupling_iteration_converged =
                        coupling_iteration_converged &&
                        _global_coupling_conv_crit[process_id]->isSatisfied();
                }
            }
            MathLib::LinAlg::copy(x, x_old);

            ++process_id;
        }  // end of for (auto& process_data : _per_process_data)

        if (coupling_iteration_converged && global_coupling_iteration > 0)
        {
            break;
        }

        if (!nonlinear_solver_status.error_norms_met)
        {
            return nonlinear_solver_status;
        }
    }

    if (!coupling_iteration_converged)
    {
        WARN(
            "The coupling iterations reaches its maximum number in time step"
            "#%u at t = %g s",
            timestep_id, t);
    }

    if (_chemical_system != nullptr)
    {
        // Sequential non-iterative approach applied here to perform water
        // chemistry calculation followed by resolving component transport
        // process.
        // TODO: move into a global loop to consider both mass balance over
        // space and localized chemical equilibrium between solutes.
        _chemical_system->doWaterChemistryCalculation(_process_solutions, dt);
    }

    int process_id = 0;
    for (auto& process_data : _per_process_data)
    {
        if (process_data->skip_time_stepping)
        {
            ++process_id;
            continue;
        }
        CoupledSolutionsForStaggeredScheme coupled_solutions(
            _solutions_of_coupled_processes, dt, process_id);

        process_data->process.setCoupledSolutionsForStaggeredScheme(
            &coupled_solutions);

        auto& pcs = process_data->process;
        auto& x = *_process_solutions[process_id];
        pcs.postTimestep(x, t, dt, process_id);
        pcs.computeSecondaryVariable(t, x, process_id);

        ++process_id;
    }

    return nonlinear_solver_status;
}

template <typename OutputClass, typename OutputClassMember>
void TimeLoop::outputSolutions(bool const output_initial_condition,
                               bool const is_staggered_coupling,
                               unsigned timestep, const double t,
                               OutputClass& output_object,
                               OutputClassMember output_class_member) const
{
    unsigned process_id = 0;
    for (auto& process_data : _per_process_data)
    {
        auto& pcs = process_data->process;
        // If nonlinear solver diverged, the solution has already been
        // saved.
        if ((!process_data->nonlinear_solver_status.error_norms_met) ||
            process_data->skip_time_stepping)
        {
            ++process_id;
            continue;
        }

        auto const& x = *_process_solutions[process_id];

        if (output_initial_condition)
        {
            pcs.preTimestep(x, _start_time,
                            process_data->timestepper->getTimeStep().dt(),
                            process_id);
            // Update secondary variables, which might be uninitialized, before
            // output.
            pcs.computeSecondaryVariable(_start_time, x, process_id);
        }
        if (is_staggered_coupling)
        {
            CoupledSolutionsForStaggeredScheme coupled_solutions(
                _solutions_of_coupled_processes, 0.0, process_id);

            process_data->process.setCoupledSolutionsForStaggeredScheme(
                &coupled_solutions);
            process_data->process
                .setCoupledTermForTheStaggeredSchemeToLocalAssemblers();
            (output_object.*output_class_member)(pcs, process_id, timestep, t,
                                                 x);
        }
        else
        {
            (output_object.*output_class_member)(pcs, process_id, timestep, t,
                                                 x);
        }

        ++process_id;
    }
}

TimeLoop::~TimeLoop()
{
    for (auto* x : _process_solutions)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }

    for (auto* x : _solutions_of_last_cpl_iteration)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }
}

}  // namespace ProcessLib
