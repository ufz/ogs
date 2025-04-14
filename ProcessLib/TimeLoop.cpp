/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TimeLoop.h"

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/contains.hpp>

#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "NumLib/ODESolver/PETScNonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "NumLib/StaggeredCoupling/StaggeredCoupling.h"
#include "NumLib/TimeStepping/Time.h"
#include "ProcessData.h"

namespace
{
void updateDeactivatedSubdomains(
    std::vector<std::unique_ptr<ProcessLib::ProcessData>> const&
        per_process_data,
    double const t)
{
    for (auto& process_data : per_process_data)
    {
        process_data->process.updateDeactivatedSubdomains(
            t, process_data->process_id);
    }
}

bool isOutputStep(std::vector<ProcessLib::Output> const& outputs,
                  const int timestep, const NumLib::Time& t,
                  const NumLib::Time& end_time)
{
    if (end_time == t)
    {
        // the last timestep is an output step
        return true;
    }

    return ranges::any_of(outputs, [timestep, t](auto const& output)
                          { return output.isOutputStep(timestep, t); });
}

void preOutputForAllProcesses(
    int const timestep, NumLib::Time const& t, double const dt,
    const NumLib::Time& end_time,
    std::vector<std::unique_ptr<ProcessLib::ProcessData>> const&
        per_process_data,
    std::vector<GlobalVector*> const& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev,
    std::vector<ProcessLib::Output> const& outputs)
{
    if (!isOutputStep(outputs, timestep, t, end_time))
    {
        return;
    }

    for (auto& process_data : per_process_data)
    {
        auto const process_id = process_data->process_id;
        auto& pcs = process_data->process;

        pcs.preOutput(t(), dt, process_solutions, process_solutions_prev,
                      process_id);
    }
}
}  // namespace

namespace ProcessLib
{
void preTimestepForAllProcesses(
    NumLib::Time const& t, double const dt,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<GlobalVector*> const& _process_solutions)
{
    for (auto& process_data : per_process_data)
    {
        auto const process_id = process_data->process_id;
        auto& pcs = process_data->process;
        pcs.preTimestep(_process_solutions, t(), dt, process_id);
    }
}

void postTimestepForAllProcesses(
    NumLib::Time const& t, double const dt,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<GlobalVector*> const& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev)
{
    for (auto& process_data : per_process_data)
    {
        auto const process_id = process_data->process_id;
        auto& pcs = process_data->process;

        pcs.computeSecondaryVariable(t(), dt, process_solutions,
                                     *process_solutions_prev[process_id],
                                     process_id);
        pcs.postTimestep(process_solutions, process_solutions_prev, t(), dt,
                         process_id);
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
    NumLib::Time const& t0,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data)
{
    std::vector<GlobalVector*> process_solutions;
    std::vector<GlobalVector*> process_solutions_prev;

    for (auto const& process_data : per_process_data)
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

    for (auto const& process_data : per_process_data)
    {
        auto& pcs = process_data->process;
        auto const process_id = process_data->process_id;
        pcs.setInitialConditions(process_solutions, process_solutions_prev,
                                 t0(), process_id);

        auto& time_disc = *process_data->time_disc;
        time_disc.setInitialState(t0());  // push IC
    }

    return {process_solutions, process_solutions_prev};
}

void calculateNonEquilibriumInitialResiduum(
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<GlobalVector*> const& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev)
{
    for (auto const& process_data : per_process_data)
    {
        auto& nonlinear_solver = process_data->nonlinear_solver;

        setEquationSystem(*process_data);
        nonlinear_solver.calculateNonEquilibriumInitialResiduum(
            process_solutions, process_solutions_prev,
            process_data->process_id);
    }
}

NumLib::NonlinearSolverStatus solveOneTimeStepOneProcess(
    std::vector<GlobalVector*>& x, std::vector<GlobalVector*> const& x_prev,
    std::size_t const timestep, double const t, double const delta_t,
    ProcessData const& process_data, std::vector<Output> const& outputs)
{
    auto& process = process_data.process;
    int const process_id = process_data.process_id;
    auto& time_disc = *process_data.time_disc;
    auto& nonlinear_solver = process_data.nonlinear_solver;

    setEquationSystem(process_data);

    // Note: Order matters!
    // First advance to the next timestep, then set known solutions at that
    // time, afterwards pass the right solution vector and time to the
    // preTimestep() hook.

    time_disc.nextTimestep(t, delta_t);

    auto const post_iteration_callback =
        [&](int const iteration, bool const converged,
            std::vector<GlobalVector*> const& x)
    {
        // Note: We don't call the postNonLinearSolver(), preOutput(),
        // computeSecondaryVariable() and postTimestep() hooks here. This might
        // lead to some inconsistencies in the data compared to regular output.
        for (auto const& output : outputs)
        {
            output.doOutputNonlinearIteration(process, process_id, timestep,
                                              NumLib::Time(t), iteration,
                                              converged, x);
        }
    };

    auto const nonlinear_solver_status =
        nonlinear_solver.solve(x, x_prev, post_iteration_callback, process_id);

    if (!nonlinear_solver_status.error_norms_met)
    {
        return nonlinear_solver_status;
    }

    process.postNonLinearSolver(x, x_prev, t, delta_t, process_id);

    return nonlinear_solver_status;
}

TimeLoop::TimeLoop(
    std::vector<Output>&& outputs,
    std::vector<std::unique_ptr<ProcessData>>&& per_process_data,
    std::unique_ptr<NumLib::StaggeredCoupling>&& staggered_coupling,
    const NumLib::Time& start_time, const NumLib::Time& end_time)
    : _outputs{std::move(outputs)},
      _per_process_data(std::move(per_process_data)),
      _start_time(start_time),
      _end_time(end_time),
      _staggered_coupling(std::move(staggered_coupling))
{
}

bool computationOfChangeNeeded(
    NumLib::TimeStepAlgorithm const& timestep_algorithm,
    NumLib::Time const& time)
{
    // for the first time step we can't compute the changes to the previous
    // time step
    if (time == timestep_algorithm.begin())
    {
        return false;
    }
    return timestep_algorithm.isSolutionErrorComputationNeeded();
}

std::pair<NumLib::TimeIncrement, bool> TimeLoop::computeTimeStepping(
    const double prev_dt, NumLib::Time& t, std::size_t& accepted_steps,
    std::size_t& rejected_steps,
    std::vector<TimeStepConstraintCallback> const& time_step_constraints)
{
    bool all_process_steps_accepted = true;
    // Get minimum time step size among step sizes of all processes.
    NumLib::TimeIncrement dt{std::numeric_limits<double>::max()};
    constexpr double eps = std::numeric_limits<double>::epsilon();

    bool const is_initial_step =
        std::any_of(_per_process_data.begin(), _per_process_data.end(),
                    [](auto const& ppd) -> bool
                    { return ppd->timestep_current.timeStepNumber() == 0; });

    for (std::size_t i = 0; i < _per_process_data.size(); i++)
    {
        auto& ppd = *_per_process_data[i];
        auto& timestep_algorithm = *ppd.timestep_algorithm.get();

        auto const& x = *_process_solutions[i];
        auto const& x_prev = *_process_solutions_prev[i];

        const double solution_error =
            computationOfChangeNeeded(timestep_algorithm, t)
                ? MathLib::LinAlg::computeRelativeNorm(
                      x, x_prev,
                      ppd.conv_crit.get() ? ppd.conv_crit->getVectorNormType()
                                          : MathLib::VecNormType::NORM2)
                : 0.0;

        ppd.timestep_current.setAccepted(
            ppd.nonlinear_solver_status.error_norms_met);

        auto [previous_step_accepted, timestepper_dt] = timestep_algorithm.next(
            solution_error, ppd.nonlinear_solver_status.number_iterations,
            ppd.timestep_previous, ppd.timestep_current);

        if (!previous_step_accepted)
        {
            // Not all processes have accepted steps.
            all_process_steps_accepted = false;
        }

        if (!ppd.nonlinear_solver_status.error_norms_met)
        {
            WARN(
                "Time step will be rejected due to nonlinear solver "
                "divergence.");
            all_process_steps_accepted = false;
        }

        if (timestepper_dt > eps || t < timestep_algorithm.end())
        {
            dt = NumLib::TimeIncrement{std::min(timestepper_dt, dt())};
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

    bool last_step_rejected = false;
    if (!is_initial_step)
    {
        if (all_process_steps_accepted)
        {
            accepted_steps++;
            last_step_rejected = false;
        }
        else
        {
            if (t <= _end_time)
            {
                t -= prev_dt;
                rejected_steps++;
                last_step_rejected = true;
            }
        }
    }

    // adjust step size considering external communciation_point_calculators
    for (auto const& time_step_constraint : time_step_constraints)
    {
        dt = NumLib::TimeIncrement{
            std::min(dt(), time_step_constraint(t, dt()))};
    }

    // Check whether the time stepping is stabilized
    if (std::abs(dt() - prev_dt) < eps)
    {
        if (last_step_rejected)
        {
            OGS_FATAL(
                "The new step size of {} is the same as that of the previous "
                "rejected time step. \nPlease re-run ogs with a proper "
                "adjustment in the numerical settings, \ne.g. those for time "
                "stepper, local or global non-linear solver.",
                dt);
        }
        else
        {
            DBUG("The time stepping is stabilized with the step size of {}.",
                 dt);
        }
    }

    // Reset the time step with the minimum step size, dt
    // Update the solution of the previous time step.
    for (std::size_t i = 0; i < _per_process_data.size(); i++)
    {
        if (all_process_steps_accepted)
        {
            auto& ppd = *_per_process_data[i];
            NumLib::updateTimeSteps(dt(), ppd.timestep_previous,
                                    ppd.timestep_current);
            auto& timestep_algorithm = ppd.timestep_algorithm;
            timestep_algorithm->resetCurrentTimeStep(
                dt(), ppd.timestep_previous, ppd.timestep_current);
        }

        auto& x = *_process_solutions[i];
        auto& x_prev = *_process_solutions_prev[i];
        if (all_process_steps_accepted)
        {
            MathLib::LinAlg::copy(x, x_prev);  // pushState
        }
        else
        {
            if (t <= _end_time)
            {
                WARN(
                    "Time step {:d} was rejected {:d} times and it will be "
                    "repeated with a reduced step size.",
                    accepted_steps + 1, _repeating_times_of_rejected_step);
                MathLib::LinAlg::copy(x_prev, x);  // popState
            }
        }
    }

    return {dt, last_step_rejected};
}

std::vector<TimeLoop::TimeStepConstraintCallback>
TimeLoop::generateOutputTimeStepConstraints(
    std::vector<double>&& fixed_times) const
{
    std::vector<TimeStepConstraintCallback> const time_step_constraints{
        [fixed_times = std::move(fixed_times)](NumLib::Time const& t, double dt)
        { return NumLib::possiblyClampDtToNextFixedTime(t, dt, fixed_times); },
        [this](NumLib::Time const& t, double dt) -> double
        {
            if (t < _end_time && _end_time < t + dt)
            {
                return _end_time() - t();
            }
            return dt;
        }};
    return time_step_constraints;
}

/// initialize output, convergence criterion, etc.
void TimeLoop::initialize()
{
    for (auto const& process_data : _per_process_data)
    {
        auto& pcs = process_data->process;
        for (auto& output : _outputs)
        {
            output.addProcess(pcs);
        }

        setTimeDiscretizedODESystem(*process_data);

        if (auto* conv_crit =
                dynamic_cast<NumLib::ConvergenceCriterionPerComponent*>(
                    process_data->conv_crit.get()))
        {
            int const process_id = process_data->process_id;
            conv_crit->setDOFTable(pcs.getDOFTable(process_id), pcs.getMesh());
        }
    }

    // initial solution storage
    std::tie(_process_solutions, _process_solutions_prev) =
        setInitialConditions(_start_time, _per_process_data);

    if (_staggered_coupling)
    {
        _staggered_coupling->initializeCoupledSolutions(_process_solutions);
    }

    updateDeactivatedSubdomains(_per_process_data, _start_time());

    auto const time_step_constraints = generateOutputTimeStepConstraints(
        calculateUniqueFixedTimesForAllOutputs(_outputs));

    std::tie(_dt, _last_step_rejected) =
        computeTimeStepping(0.0, _current_time, _accepted_steps,
                            _rejected_steps, time_step_constraints);

    // Output initial conditions
    {
        preOutputInitialConditions(_start_time, _dt());
        outputSolutions(0, _start_time(), &Output::doOutput);
    }

    calculateNonEquilibriumInitialResiduum(
        _per_process_data, _process_solutions, _process_solutions_prev);
}

bool TimeLoop::executeTimeStep()
{
    BaseLib::RunTime time_timestep;
    time_timestep.start();

    _current_time += _dt();

    const std::size_t timesteps = _accepted_steps + 1;
    // TODO(wenqing): , input option for time unit.
    INFO(
        "Time stepping at step #{:d} and time {} with step size {} "
        "started.",
        timesteps, _current_time, _dt);

    updateDeactivatedSubdomains(_per_process_data, _current_time());

    successful_time_step =
        preTsNonlinearSolvePostTs(_current_time, _dt(), timesteps);
    INFO("[time] Time step #{:d} took {:g} s.", timesteps,
         time_timestep.elapsed());
    return successful_time_step;
}

bool TimeLoop::calculateNextTimeStep()
{
    const double prev_dt = _dt();
    auto const current_time = _current_time;

    const std::size_t timesteps = _accepted_steps + 1;

    auto const time_step_constraints = generateOutputTimeStepConstraints(
        calculateUniqueFixedTimesForAllOutputs(_outputs));

    // _last_step_rejected is also checked in computeTimeStepping.
    std::tie(_dt, _last_step_rejected) =
        computeTimeStepping(prev_dt, _current_time, _accepted_steps,
                            _rejected_steps, time_step_constraints);

    if (!_last_step_rejected)
    {
        outputSolutions(timesteps, current_time(), &Output::doOutput);
    }

    if (current_time == (_current_time + _dt()))
    {
        ERR("Time step size of {} is too small.\n"
            "Time stepping stops at step {:d} and at time of {}.",
            _dt, timesteps, _current_time);
        return false;
    }

    if (_current_time >= _end_time || _current_time + _dt() > _end_time)
    {
        return false;
    }

    return true;
}

void TimeLoop::outputLastTimeStep() const
{
    INFO(
        "The whole computation of the time stepping took {:d} steps, in which\n"
        "\t the accepted steps are {:d}, and the rejected steps are {:d}.\n",
        _accepted_steps + _rejected_steps, _accepted_steps, _rejected_steps);

    // output last time step
    if (successful_time_step)
    {
        outputSolutions(_accepted_steps + _rejected_steps, _current_time(),
                        &Output::doOutputLastTimestep);
    }
}

bool TimeLoop::preTsNonlinearSolvePostTs(NumLib::Time const& t, double const dt,
                                         std::size_t const timesteps)
{
    preTimestepForAllProcesses(t, dt, _per_process_data, _process_solutions);

    NumLib::NonlinearSolverStatus nonlinear_solver_status;

    if (_staggered_coupling)
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
        // Later on, the timestep_algorithm might reject the timestep. We assume
        // that this is a rare case, so still, we call preOutput() here. We
        // don't expect a large overhead from it.
        preOutputForAllProcesses(timesteps, t, dt, _end_time, _per_process_data,
                                 _process_solutions, _process_solutions_prev,
                                 _outputs);

        postTimestepForAllProcesses(t, dt, _per_process_data,
                                    _process_solutions,
                                    _process_solutions_prev);
    }
    return nonlinear_solver_status.error_norms_met;
}

static NumLib::NonlinearSolverStatus solveMonolithicProcess(
    const NumLib::Time& t, const double dt, const std::size_t timestep_id,
    ProcessData const& process_data, std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev,
    std::vector<Output> const& outputs)
{
    BaseLib::RunTime time_timestep_process;
    time_timestep_process.start();

    //    INFO(
    //        "Solving process #{:d} at time step #{:d} and time {:f} with step
    //        size {:f} started.", process_data.process_id, timestep_id, t(),
    //        dt);
    //

    INFO("Solving process #{:d} started.", process_data.process_id);

    auto const nonlinear_solver_status = solveOneTimeStepOneProcess(
        x, x_prev, timestep_id, t(), dt, process_data, outputs);

    INFO("[time] Solving process #{:d} took {:g} s in time step #{:d}",
         process_data.process_id, time_timestep_process.elapsed(), timestep_id);

    return nonlinear_solver_status;
}

static constexpr std::string_view timestepper_cannot_reduce_dt =
    "Time stepper cannot reduce the time step size further.";

NumLib::NonlinearSolverStatus TimeLoop::solveUncoupledEquationSystems(
    const NumLib::Time& t, const double dt, const std::size_t timestep_id)
{
    NumLib::NonlinearSolverStatus nonlinear_solver_status;

    for (auto const& process_data : _per_process_data)
    {
        auto const process_id = process_data->process_id;
        nonlinear_solver_status = solveMonolithicProcess(
            t, dt, timestep_id, *process_data, _process_solutions,
            _process_solutions_prev, _outputs);

        process_data->nonlinear_solver_status = nonlinear_solver_status;
        if (!nonlinear_solver_status.error_norms_met)
        {
            ERR("The nonlinear solver failed in time step #{:d} at t = {} s "
                "for process #{:d}.",
                timestep_id, t, process_id);

            if (!process_data->timestep_algorithm->canReduceTimestepSize(
                    process_data->timestep_current,
                    process_data->timestep_previous))
            {
                // save unsuccessful solution
                for (auto const& output : _outputs)
                {
                    output.doOutputAlways(
                        process_data->process, process_id, timestep_id, t,
                        process_data->nonlinear_solver_status.number_iterations,
                        process_data->nonlinear_solver_status.error_norms_met,
                        _process_solutions);
                }
                OGS_FATAL(timestepper_cannot_reduce_dt.data());
            }

            return nonlinear_solver_status;
        }
    }

    return nonlinear_solver_status;
}

NumLib::NonlinearSolverStatus
TimeLoop::solveCoupledEquationSystemsByStaggeredScheme(
    const NumLib::Time& t, const double dt, const std::size_t timestep_id)
{
    auto const nonlinear_solver_status =
        _staggered_coupling->execute<ProcessData, Output>(
            t(), dt, timestep_id, _process_solutions, _process_solutions_prev,
            _per_process_data, _outputs, &solveOneTimeStepOneProcess);

    _last_step_rejected = nonlinear_solver_status.error_norms_met;

    {
        for (auto const& process_data : _per_process_data)
        {
            auto& pcs = process_data->process;
            int const process_id = process_data->process_id;
            auto& ode_sys = *process_data->tdisc_ode_sys;
            pcs.solveReactionEquation(_process_solutions,
                                      _process_solutions_prev, t(), dt, ode_sys,
                                      process_id);
        }
    }

    return nonlinear_solver_status;
}

template <typename OutputClassMember>
void TimeLoop::outputSolutions(unsigned timestep, const double t,
                               OutputClassMember output_class_member) const
{
    for (auto const& process_data : _per_process_data)
    {
        // If nonlinear solver diverged, the solution has already been
        // saved.
        if (!process_data->nonlinear_solver_status.error_norms_met)
        {
            continue;
        }

        auto const process_id = process_data->process_id;
        auto const& pcs = process_data->process;

        for (auto const& output_object : _outputs)
        {
            (output_object.*output_class_member)(
                pcs, process_id, timestep, NumLib::Time(t),
                process_data->nonlinear_solver_status.number_iterations,
                process_data->nonlinear_solver_status.error_norms_met,
                _process_solutions);
        }
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
}

void TimeLoop::preOutputInitialConditions(NumLib::Time const& t,
                                          const double dt) const
{
    for (auto const& process_data : _per_process_data)
    {
        // If nonlinear solver diverged, the solution has already been
        // saved.
        if (!process_data->nonlinear_solver_status.error_norms_met)
        {
            continue;
        }

        auto const process_id = process_data->process_id;
        auto& pcs = process_data->process;

        process_data->time_disc->nextTimestep(t(), dt);

        pcs.preTimestep(_process_solutions, _start_time(), dt, process_id);

        pcs.preOutput(_start_time(), dt, _process_solutions,
                      _process_solutions_prev, process_id);

        // Update secondary variables, which might be uninitialized, before
        // output.
        pcs.computeSecondaryVariable(_start_time(), dt, _process_solutions,
                                     *_process_solutions_prev[process_id],
                                     process_id);
    }
}
}  // namespace ProcessLib
