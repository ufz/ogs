/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TimeLoop.h"

#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "CoupledSolutionsForStaggeredScheme.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "NumLib/ODESolver/PETScNonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "ProcessData.h"
#include "ProcessLib/CreateProcessData.h"
#include "ProcessLib/Output/CreateOutput.h"
namespace
{
void setEquationSystem(NumLib::NonlinearSolverBase& nonlinear_solver,
                       NumLib::EquationSystem& eq_sys,
                       NumLib::ConvergenceCriterion& conv_crit,
                       NumLib::NonlinearSolverTag nl_tag)
{
    using Tag = NumLib::NonlinearSolverTag;
    switch (nl_tag)
    {
        case Tag::Picard:
        {
            using EqSys = NumLib::NonlinearSystem<Tag::Picard>;
            auto& eq_sys_ = static_cast<EqSys&>(eq_sys);
            if (auto* nl_solver =
                    dynamic_cast<NumLib::NonlinearSolver<Tag::Picard>*>(
                        &nonlinear_solver);
                nl_solver != nullptr)
            {
                nl_solver->setEquationSystem(eq_sys_, conv_crit);
            }
            else
            {
                OGS_FATAL(
                    "Could not cast nonlinear solver to Picard type solver.");
            }
            break;
        }
        case Tag::Newton:
        {
            using EqSys = NumLib::NonlinearSystem<Tag::Newton>;
            auto& eq_sys_ = static_cast<EqSys&>(eq_sys);

            if (auto* nl_solver =
                    dynamic_cast<NumLib::NonlinearSolver<Tag::Newton>*>(
                        &nonlinear_solver);
                nl_solver != nullptr)
            {
                nl_solver->setEquationSystem(eq_sys_, conv_crit);
            }
#ifdef USE_PETSC
            else if (auto* nl_solver =
                         dynamic_cast<NumLib::PETScNonlinearSolver*>(
                             &nonlinear_solver);
                     nl_solver != nullptr)
            {
                nl_solver->setEquationSystem(eq_sys_, conv_crit);
            }
#endif  // USE_PETSC
            else
            {
                OGS_FATAL(
                    "Could not cast nonlinear solver to Newton type solver.");
            }
            break;
        }
    }
}

bool isMonolithicProcess(ProcessLib::ProcessData const& process_data)
{
    return process_data.process.isMonolithicSchemeUsed();
}
}  // namespace

namespace ProcessLib
{
void preTimestepForAllProcesses(
    double const t, double const dt,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<GlobalVector*> const& _process_solutions)
{
    for (auto& process_data : per_process_data)
    {
        auto const process_id = process_data->process_id;
        auto& pcs = process_data->process;
        pcs.preTimestep(_process_solutions, t, dt, process_id);
    }
}

void postTimestepForAllProcesses(
    double const t, double const dt,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<GlobalVector*> const& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev)
{
    std::vector<GlobalVector*> x_dots;
    x_dots.reserve(per_process_data.size());
    for (auto& process_data : per_process_data)
    {
        auto const process_id = process_data->process_id;
        auto const& ode_sys = *process_data->tdisc_ode_sys;
        auto const& time_discretization = *process_data->time_disc;

        x_dots.emplace_back(&NumLib::GlobalVectorProvider::provider.getVector(
            ode_sys.getMatrixSpecifications(process_id)));

        time_discretization.getXdot(*process_solutions[process_id],
                                    *process_solutions_prev[process_id],
                                    *x_dots[process_id]);
    }

    // All _per_process_data share the first process.
    bool const is_staggered_coupling =
        !isMonolithicProcess(*per_process_data[0]);

    for (auto& process_data : per_process_data)
    {
        auto const process_id = process_data->process_id;
        auto& pcs = process_data->process;

        if (is_staggered_coupling)
        {
            CoupledSolutionsForStaggeredScheme coupled_solutions(
                process_solutions);
            pcs.setCoupledSolutionsForStaggeredScheme(&coupled_solutions);
        }
        auto& x_dot = *x_dots[process_id];
        pcs.computeSecondaryVariable(t, dt, process_solutions, x_dot,
                                     process_id);
        pcs.postTimestep(process_solutions, t, dt, process_id);
    }
    for (auto& x_dot : x_dots)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x_dot);
    }
}

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
    // TODO (naumov) Provide a function to nonlinear_solver to distinguish the
    // types. Could be handy, because a nonlinear solver could handle both types
    // like PETScSNES.
    else if ((dynamic_cast<NonlinearSolverNewton*>(
                  &process_data.nonlinear_solver) != nullptr)
#ifdef USE_PETSC
             || (dynamic_cast<NumLib::PETScNonlinearSolver*>(
                     &process_data.nonlinear_solver) != nullptr)
#endif  // USE_PETSC
    )
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
}

void setTimeDiscretizedODESystem(ProcessData& process_data)
{
    setTimeDiscretizedODESystem(process_data, process_data.process);
}

std::pair<std::vector<GlobalVector*>, std::vector<GlobalVector*>>
setInitialConditions(
    double const t0,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data)
{
    std::vector<GlobalVector*> process_solutions;
    std::vector<GlobalVector*> process_solutions_prev;

    for (auto& process_data : per_process_data)
    {
        auto const process_id = process_data->process_id;
        auto& ode_sys = *process_data->tdisc_ode_sys;

        // append a solution vector of suitable size
        process_solutions.emplace_back(
            &NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications(process_id)));
        process_solutions_prev.emplace_back(
            &NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications(process_id)));
    }

    for (auto& process_data : per_process_data)
    {
        auto& pcs = process_data->process;
        auto const process_id = process_data->process_id;
        pcs.setInitialConditions(process_solutions, process_solutions_prev, t0,
                                 process_id);

        auto& time_disc = *process_data->time_disc;
        time_disc.setInitialState(t0);     // push IC
    }

    return {process_solutions, process_solutions_prev};
}

void calculateNonEquilibriumInitialResiduum(
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<GlobalVector*>
        process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev)
{
    INFO("Calculate non-equilibrium initial residuum.");
    for (auto& process_data : per_process_data)
    {
        auto& ode_sys = *process_data->tdisc_ode_sys;

        auto& time_disc = *process_data->time_disc;
        auto& nonlinear_solver = process_data->nonlinear_solver;
        auto& conv_crit = *process_data->conv_crit;

        auto const nl_tag = process_data->nonlinear_solver_tag;
        setEquationSystem(nonlinear_solver, ode_sys, conv_crit, nl_tag);
        // dummy values to handle the time derivative terms more or less
        // correctly, i.e. to ignore them.
        double const t = 0;
        double const dt = 1;
        time_disc.nextTimestep(t, dt);
        nonlinear_solver.calculateNonEquilibriumInitialResiduum(
            process_solutions, process_solutions_prev,
            process_data->process_id);
    }
}

NumLib::NonlinearSolverStatus solveOneTimeStepOneProcess(
    std::vector<GlobalVector*>& x, std::vector<GlobalVector*> const& x_prev,
    std::size_t const timestep, double const t, double const delta_t,
    ProcessData const& process_data, Output& output_control)
{
    auto& process = process_data.process;
    int const process_id = process_data.process_id;
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

    auto const post_iteration_callback =
        [&](int iteration, std::vector<GlobalVector*> const& x) {
            output_control.doOutputNonlinearIteration(
                process, process_id, timestep, t, iteration, x);
        };

    auto const nonlinear_solver_status =
        nonlinear_solver.solve(x, x_prev, post_iteration_callback, process_id);

    if (nonlinear_solver_status.error_norms_met)
    {
        GlobalVector& x_dot = NumLib::GlobalVectorProvider::provider.getVector(
            ode_sys.getMatrixSpecifications(process_id));

        time_disc.getXdot(*x[process_id], *x_prev[process_id], x_dot);

        process.postNonLinearSolver(*x[process_id], x_dot, t, delta_t,
                                    process_id);
        NumLib::GlobalVectorProvider::provider.releaseVector(x_dot);
    }

    return nonlinear_solver_status;
}

TimeLoop::TimeLoop(std::unique_ptr<Output>&& output,
                   std::vector<std::unique_ptr<ProcessData>>&& per_process_data,
                   const int global_coupling_max_iterations,
                   std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>&&
                       global_coupling_conv_crit,
                   const double start_time, const double end_time)
    : _output(std::move(output)),
      _per_process_data(std::move(per_process_data)),
      _start_time(start_time),
      _end_time(end_time),
      _global_coupling_max_iterations(global_coupling_max_iterations),
      _global_coupling_conv_crit(std::move(global_coupling_conv_crit))
{
}

void TimeLoop::setCoupledSolutions()
{
    for (auto& process_data : _per_process_data)
    {
        auto const& x = *_process_solutions[process_data->process_id];

        // Create a vector to store the solution of the last coupling iteration
        auto& x0 = NumLib::GlobalVectorProvider::provider.getVector(x);
        MathLib::LinAlg::copy(x, x0);

        // append a solution vector of suitable size
        _solutions_of_last_cpl_iteration.emplace_back(&x0);
    }
}

double TimeLoop::computeTimeStepping(const double prev_dt, double& t,
                                     std::size_t& accepted_steps,
                                     std::size_t& rejected_steps)
{
    bool all_process_steps_accepted = true;
    // Get minimum time step size among step sizes of all processes.
    double dt = std::numeric_limits<double>::max();

    bool is_initial_step = false;
    for (std::size_t i = 0; i < _per_process_data.size(); i++)
    {
        auto& ppd = *_per_process_data[i];
        const auto& timestepper = ppd.timestepper;
        if (timestepper->getTimeStep().steps() == 0)
        {
            is_initial_step = true;
        }

        auto& time_disc = ppd.time_disc;
        auto const& x = *_process_solutions[i];
        auto const& x_prev = *_process_solutions_prev[i];

        auto const& conv_crit = ppd.conv_crit;
        const MathLib::VecNormType norm_type =
            (conv_crit) ? conv_crit->getVectorNormType()
                        : MathLib::VecNormType::NORM2;

        const double solution_error =
            (timestepper->isSolutionErrorComputationNeeded())
                ? ((t == timestepper->begin())
                       ? 0.  // Always accepts the zeroth step
                       : time_disc->computeRelativeChangeFromPreviousTimestep(
                             x, x_prev, norm_type))
                : 0.;

        if (!ppd.nonlinear_solver_status.error_norms_met)
        {
            timestepper->setAcceptedOrNot(false);
        }
        else
        {
            timestepper->setAcceptedOrNot(true);
        }

        auto [step_accepted, timestepper_dt] = timestepper->next(
            solution_error, ppd.nonlinear_solver_status.number_iterations);

        if (!step_accepted &&
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

        if (timestepper_dt > std::numeric_limits<double>::epsilon() ||
            std::abs(t - timestepper->end()) <
                std::numeric_limits<double>::epsilon())
        {
            dt = std::min(timestepper_dt, dt);
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

    dt = NumLib::possiblyClampDtToNextFixedTime(t, dt,
                                                _output->getFixedOutputTimes());
    // Check whether the time stepping is stabilized
    if (std::abs(dt - prev_dt) < std::numeric_limits<double>::epsilon())
    {
        if (_last_step_rejected)
        {
            OGS_FATAL(
                "The new step size of {:g} is the same as that of the previous "
                "rejected time step. \nPlease re-run ogs with a proper "
                "adjustment in the numerical settings, \ne.g those for "
                "time stepper, local or global non-linear solver.",
                dt);
        }
        else
        {
            DBUG("The time stepping is stabilized with the step size of {:g}.",
                 dt);
        }
    }

    // Reset the time step with the minimum step size, dt
    // Update the solution of the previous time step.
    for (std::size_t i = 0; i < _per_process_data.size(); i++)
    {
        const auto& ppd = *_per_process_data[i];
        auto& timestepper = ppd.timestepper;
        if (all_process_steps_accepted)
        {
            timestepper->resetCurrentTimeStep(dt);
        }

        if (t == timestepper->begin())
        {
            continue;
        }

        auto& x = *_process_solutions[i];
        auto& x_prev = *_process_solutions_prev[i];
        if (all_process_steps_accepted)
        {
            MathLib::LinAlg::copy(x, x_prev);  // pushState
        }
        else
        {
            if (t < _end_time || std::abs(t - _end_time) <
                                     std::numeric_limits<double>::epsilon())
            {
                WARN(
                    "Time step {:d} was rejected {:d} times "
                    "and it will be repeated with a reduced step size.",
                    accepted_steps + 1, _repeating_times_of_rejected_step);
                MathLib::LinAlg::copy(x_prev, x);  // popState
            }
        }
    }

    return dt;
}

/// initialize output, convergence criterion, etc.
void TimeLoop::initialize()
{
    for (auto& process_data : _per_process_data)
    {
        auto& pcs = process_data->process;
        int const process_id = process_data->process_id;
        _output->addProcess(pcs);

        setTimeDiscretizedODESystem(*process_data);

        if (auto* conv_crit =
                dynamic_cast<NumLib::ConvergenceCriterionPerComponent*>(
                    process_data->conv_crit.get()))
        {
            conv_crit->setDOFTable(pcs.getDOFTable(process_id), pcs.getMesh());
        }
    }

    // initial solution storage
    std::tie(_process_solutions, _process_solutions_prev) =
        setInitialConditions(_start_time, _per_process_data);

    // All _per_process_data share the first process.
    bool const is_staggered_coupling =
        !isMonolithicProcess(*_per_process_data[0]);

    if (is_staggered_coupling)
    {
        setCoupledSolutions();
    }

    // Output initial conditions
    {
        const bool output_initial_condition = true;
        outputSolutions(output_initial_condition, 0, _start_time, *_output,
                        &Output::doOutput);
    }
}

/*
 * TODO:
 * Now we have a structure inside the time loop which is very similar to the
 * nonlinear solver. And admittedly, the control flow inside the nonlinear
 * solver is rather complicated. Maybe in the future one can introduce an
 * abstraction that can do both the convergence checks of the coupling loop and
 * of the nonlinear solver.
 */
bool TimeLoop::loop()
{
    // All _per_process_data share the first process.
    bool const is_staggered_coupling =
        !isMonolithicProcess(*_per_process_data[0]);

    bool non_equilibrium_initial_residuum_computed = false;
    double t = _start_time;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    NumLib::NonlinearSolverStatus nonlinear_solver_status;

    double dt = computeTimeStepping(0.0, t, accepted_steps, rejected_steps);

    while (t < _end_time)
    {
        BaseLib::RunTime time_timestep;
        time_timestep.start();

        t += dt;
        const double prev_dt = dt;

        const std::size_t timesteps = accepted_steps + 1;
        // TODO(wenqing): , input option for time unit.
        INFO(
            "=== Time stepping at step #{:d} and time {:g} with step size {:g}",
            timesteps, t, dt);

        // Check element deactivation:
        for (auto& process_data : _per_process_data)
        {
            process_data->process.updateDeactivatedSubdomains(
                t, process_data->process_id);
        }

        if (!non_equilibrium_initial_residuum_computed)
        {
            calculateNonEquilibriumInitialResiduum(
                _per_process_data, _process_solutions, _process_solutions_prev);
            non_equilibrium_initial_residuum_computed = true;
        }

        preTimestepForAllProcesses(t, dt, _per_process_data,
                                   _process_solutions);

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

        // Run post time step only if the last iteration was successful.
        // Otherwise it runs the risks to get the same errors as in the last
        // iteration, an exception thrown in assembly, for example.
        if (nonlinear_solver_status.error_norms_met)
        {
            postTimestepForAllProcesses(t, dt, _per_process_data,
                                        _process_solutions,
                                        _process_solutions_prev);
        }

        INFO("[time] Time step #{:d} took {:g} s.", timesteps,
             time_timestep.elapsed());


        if (!_last_step_rejected)
        {
            const bool output_initial_condition = false;
            outputSolutions(output_initial_condition, timesteps, t, *_output,
                            &Output::doOutput);
        }

        dt = computeTimeStepping(prev_dt, t, accepted_steps, rejected_steps);
        if (std::abs(t - _end_time) < std::numeric_limits<double>::epsilon() ||
            t + dt > _end_time)
        {
            break;
        }

        if (dt < std::numeric_limits<double>::epsilon())
        {
            WARN(
                "Time step size of {:g} is too small.\n"
                "Time stepping stops at step {:d} and at time of {:g}.",
                dt, timesteps, t);
            break;
        }
    }

    INFO(
        "The whole computation of the time stepping took {:d} steps, in which\n"
        "\t the accepted steps are {:d}, and the rejected steps are {:d}.\n",
        accepted_steps + rejected_steps, accepted_steps, rejected_steps);

    // output last time step
    if (nonlinear_solver_status.error_norms_met)
    {
        const bool output_initial_condition = false;
        outputSolutions(output_initial_condition,
                        accepted_steps + rejected_steps, t, *_output,
                        &Output::doOutputLastTimestep);
    }

    return nonlinear_solver_status.error_norms_met;
}

static NumLib::NonlinearSolverStatus solveMonolithicProcess(
    const double t, const double dt, const std::size_t timestep_id,
    ProcessData const& process_data, std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev, Output& output)
{
    BaseLib::RunTime time_timestep_process;
    time_timestep_process.start();

    auto const nonlinear_solver_status = solveOneTimeStepOneProcess(
        x, x_prev, timestep_id, t, dt, process_data, output);

    INFO("[time] Solving process #{:d} took {:g} s in time step #{:d} ",
         process_data.process_id, time_timestep_process.elapsed(), timestep_id);

    return nonlinear_solver_status;
}

static constexpr std::string_view timestepper_cannot_reduce_dt =
    "Time stepper cannot reduce the time step size further.";

NumLib::NonlinearSolverStatus TimeLoop::solveUncoupledEquationSystems(
    const double t, const double dt, const std::size_t timestep_id)
{
    NumLib::NonlinearSolverStatus nonlinear_solver_status;
    for (auto& process_data : _per_process_data)
    {
        auto const process_id = process_data->process_id;
        nonlinear_solver_status = solveMonolithicProcess(
            t, dt, timestep_id, *process_data, _process_solutions,
            _process_solutions_prev, *_output);

        process_data->nonlinear_solver_status = nonlinear_solver_status;
        if (!nonlinear_solver_status.error_norms_met)
        {
            ERR("The nonlinear solver failed in time step #{:d} at t = {:g} s "
                "for process #{:d}.",
                timestep_id, t, process_id);

            if (!process_data->timestepper->canReduceTimestepSize())
            {
                // save unsuccessful solution
                _output->doOutputAlways(
                    process_data->process, process_id, timestep_id, t,
                    process_data->nonlinear_solver_status.number_iterations,
                    _process_solutions);
                OGS_FATAL(timestepper_cannot_reduce_dt.data());
            }

            return nonlinear_solver_status;
        }
    }

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

    NumLib::NonlinearSolverStatus nonlinear_solver_status{false, -1};
    bool coupling_iteration_converged = true;
    for (int global_coupling_iteration = 0;
         global_coupling_iteration < _global_coupling_max_iterations;
         global_coupling_iteration++, resetCouplingConvergenceCriteria())
    {
        // TODO(wenqing): use process name
        coupling_iteration_converged = true;
        int const last_process_id = _per_process_data.size() - 1;
        for (auto& process_data : _per_process_data)
        {
            auto const process_id = process_data->process_id;
            BaseLib::RunTime time_timestep_process;
            time_timestep_process.start();

            // The following setting of coupled_solutions can be removed only if
            // the CoupledSolutionsForStaggeredScheme and related functions are
            // removed totally from the computation of the secondary variable
            // and from post-time functions.
            CoupledSolutionsForStaggeredScheme coupled_solutions(
                _process_solutions);

            process_data->process.setCoupledSolutionsForStaggeredScheme(
                &coupled_solutions);

            nonlinear_solver_status = solveOneTimeStepOneProcess(
                _process_solutions, _process_solutions_prev, timestep_id, t, dt,
                *process_data, *_output);
            process_data->nonlinear_solver_status = nonlinear_solver_status;

            INFO(
                "[time] Solving process #{:d} took {:g} s in time step #{:d} "
                " coupling iteration #{:d}",
                process_id, time_timestep_process.elapsed(), timestep_id,
                global_coupling_iteration);

            if (!nonlinear_solver_status.error_norms_met)
            {
                WARN(
                    "The nonlinear solver failed in time step #{:d} at t = "
                    "{:g} s for process #{:d}.",
                    timestep_id, t, process_id);
                _last_step_rejected = true;
                return nonlinear_solver_status;
            }

            // Check the convergence of the coupling iteration
            auto& x = *_process_solutions[process_id];
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
            "#{:d} at t = {:g} s",
            timestep_id, t);
    }

    {
        auto& pcs = _per_process_data[0]->process;
        pcs.solveReactionEquation(_process_solutions, t, dt);
    }

    return nonlinear_solver_status;
}

template <typename OutputClass, typename OutputClassMember>
void TimeLoop::outputSolutions(bool const output_initial_condition,
                               unsigned timestep, const double t,
                               OutputClass& output_object,
                               OutputClassMember output_class_member) const
{
    // All _per_process_data share the first process.
    bool const is_staggered_coupling =
        !isMonolithicProcess(*_per_process_data[0]);

    for (auto& process_data : _per_process_data)
    {
        // If nonlinear solver diverged, the solution has already been
        // saved.
        if (!process_data->nonlinear_solver_status.error_norms_met)
        {
            continue;
        }

        auto const process_id = process_data->process_id;
        auto& pcs = process_data->process;

        if (!is_staggered_coupling && output_initial_condition)
        {
            auto const& ode_sys = *process_data->tdisc_ode_sys;
            // dummy value to handle the time derivative terms more or less
            // correctly, i.e. to ignore them.
            double const dt = 1;
            process_data->time_disc->nextTimestep(t, dt);

            auto& x_dot = NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications(process_id));
            x_dot.setZero();

            pcs.preTimestep(_process_solutions, _start_time, dt, process_id);
            // Update secondary variables, which might be uninitialized, before
            // output.
            pcs.computeSecondaryVariable(_start_time, dt, _process_solutions,
                                         x_dot, process_id);

            NumLib::GlobalVectorProvider::provider.releaseVector(x_dot);
        }
        if (is_staggered_coupling && output_initial_condition)
        {
            CoupledSolutionsForStaggeredScheme coupled_solutions(
                _process_solutions);

            process_data->process.setCoupledSolutionsForStaggeredScheme(
                &coupled_solutions);
            process_data->process
                .setCoupledTermForTheStaggeredSchemeToLocalAssemblers(
                    process_id);

            auto const& ode_sys = *process_data->tdisc_ode_sys;
            // dummy value to handle the time derivative terms more or less
            // correctly, i.e. to ignore them.
            double const dt = 1;
            process_data->time_disc->nextTimestep(t, dt);

            auto& x_dot = NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications(process_id));
            x_dot.setZero();

            pcs.preTimestep(_process_solutions, _start_time, dt, process_id);
            // Update secondary variables, which might be uninitialized, before
            // output.
            pcs.computeSecondaryVariable(_start_time, dt, _process_solutions,
                                         x_dot, process_id);

            NumLib::GlobalVectorProvider::provider.releaseVector(x_dot);
        }
        (output_object.*output_class_member)(
            pcs, process_id, timestep, t,
            process_data->nonlinear_solver_status.number_iterations,
            _process_solutions);
    }
}

TimeLoop::~TimeLoop()
{
    for (auto* x : _process_solutions)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }
    for (auto* x : _process_solutions_prev)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }

    for (auto* x : _solutions_of_last_cpl_iteration)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }
}

}  // namespace ProcessLib
