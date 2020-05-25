/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TimeLoop.h"

#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "ChemistryLib/ChemicalSolverInterface.h"
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

bool isMonolithicProcess(ProcessLib::ProcessData const& process_data)
{
    return process_data.process.isMonolithicSchemeUsed();
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
        auto& pcs = process_data->process;
        auto const process_id = process_data->process_id;
        auto& time_disc = *process_data->time_disc;

        auto& ode_sys = *process_data->tdisc_ode_sys;

        // append a solution vector of suitable size
        process_solutions.emplace_back(
            &NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications(process_id)));
        process_solutions_prev.emplace_back(
            &NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications(process_id)));

        auto& x = *process_solutions[process_id];
        auto& x_prev = *process_solutions_prev[process_id];
        pcs.setInitialConditions(process_id, t0, x);
        MathLib::LinAlg::finalizeAssembly(x);

        time_disc.setInitialState(t0);     // push IC
        MathLib::LinAlg::copy(x, x_prev);  // pushState
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
                process, process_id, timestep, t, x, iteration);
        };

    auto const nonlinear_solver_status =
        nonlinear_solver.solve(x, x_prev, post_iteration_callback, process_id);

    if (nonlinear_solver_status.error_norms_met)
    {
        process.postNonLinearSolver(*x[process_id], t, delta_t, process_id);
    }

    return nonlinear_solver_status;
}

TimeLoop::TimeLoop(std::unique_ptr<Output>&& output,
                   std::vector<std::unique_ptr<ProcessData>>&& per_process_data,
                   const int global_coupling_max_iterations,
                   std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>&&
                       global_coupling_conv_crit,
                   std::unique_ptr<ChemistryLib::ChemicalSolverInterface>&&
                       chemical_solver_interface,
                   const double start_time, const double end_time)
    : output_(std::move(output)),
      per_process_data_(std::move(per_process_data)),
      start_time_(start_time),
      end_time_(end_time),
      global_coupling_max_iterations_(global_coupling_max_iterations),
      global_coupling_conv_crit_(std::move(global_coupling_conv_crit)),
      chemical_solver_interface_(std::move(chemical_solver_interface))
{
}

void TimeLoop::setCoupledSolutions()
{
    for (auto& process_data : per_process_data_)
    {
        auto const& x = *process_solutions_[process_data->process_id];

        // Create a vector to store the solution of the last coupling iteration
        auto& x0 = NumLib::GlobalVectorProvider::provider.getVector(x);
        MathLib::LinAlg::copy(x, x0);

        // append a solution vector of suitable size
        solutions_of_last_cpl_iteration_.emplace_back(&x0);
    }
}

double TimeLoop::computeTimeStepping(const double prev_dt, double& t,
                                     std::size_t& accepted_steps,
                                     std::size_t& rejected_steps)
{
    bool all_process_steps_accepted = true;
    // Get minimum time step size among step sizes of all processes.
    double dt = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < per_process_data_.size(); i++)
    {
        auto& ppd = *per_process_data_[i];
        const auto& timestepper = ppd.timestepper;

        auto& time_disc = ppd.time_disc;
        auto const& x = *process_solutions_[i];
        auto const& x_prev = *process_solutions_prev_[i];

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
    }

    if (all_process_steps_accepted)
    {
        repeating_times_of_rejected_step_ = 0;
    }
    else
    {
        repeating_times_of_rejected_step_++;
    }

    bool is_initial_step = false;
    // Reset the time step with the minimum step size, dt
    // Update the solution of the previous time step.
    for (std::size_t i = 0; i < per_process_data_.size(); i++)
    {
        const auto& ppd = *per_process_data_[i];
        auto& timestepper = ppd.timestepper;
        timestepper->resetCurrentTimeStep(dt);

        if (t == timestepper->begin())
        {
            is_initial_step = true;
            continue;
        }

        auto& x = *process_solutions_[i];
        auto& x_prev = *process_solutions_prev_[i];
        if (all_process_steps_accepted)
        {
            MathLib::LinAlg::copy(x, x_prev);  // pushState
        }
        else
        {
            if (t < end_time_ || std::abs(t - end_time_) <
                                     std::numeric_limits<double>::epsilon())
            {
                WARN(
                    "Time step {:d} was rejected {:d} times "
                    "and it will be repeated with a reduced step size.",
                    accepted_steps + 1, repeating_times_of_rejected_step_);
                MathLib::LinAlg::copy(x_prev, x);  // popState
            }
        }
    }

    if (!is_initial_step)
    {
        if (all_process_steps_accepted)
        {
            accepted_steps++;
            last_step_rejected_ = false;
        }
        else
        {
            if (t < end_time_ || std::abs(t - end_time_) <
                                     std::numeric_limits<double>::epsilon())
            {
                t -= prev_dt;
                rejected_steps++;
                last_step_rejected_ = true;
            }
        }
    }

    // Adjust step size if t < end_time_, while t+dt exceeds the end time
    if (t < end_time_ && t + dt > end_time_)
    {
        dt = end_time_ - t;
    }

    // Check whether the time stepping is stabilized
    if (std::fabs(dt - prev_dt) < std::numeric_limits<double>::epsilon())
    {
        if (last_step_rejected_)
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

    return dt;
}

/// initialize output, convergence criterion, etc.
void TimeLoop::initialize()
{
    for (auto& process_data : per_process_data_)
    {
        auto& pcs = process_data->process;
        int const process_id = process_data->process_id;
        output_->addProcess(pcs, process_id);

        setTimeDiscretizedODESystem(*process_data);

        if (auto* conv_crit =
                dynamic_cast<NumLib::ConvergenceCriterionPerComponent*>(
                    process_data->conv_crit.get()))
        {
            conv_crit->setDOFTable(pcs.getDOFTable(process_id), pcs.getMesh());
        }

        // Add the fixed times of output to time stepper in order that
        // the time stepping is performed and the results are output at
        // these times. Note: only the adaptive time steppers can have the
        // the fixed times.
        auto& timestepper = process_data->timestepper;
        timestepper->addFixedOutputTimes(output_->getFixedOutputTimes());
    }

    // init solution storage
    std::tie(process_solutions_, process_solutions_prev_) =
        setInitialConditions(start_time_, per_process_data_);

    if (chemical_solver_interface_)
    {
        BaseLib::RunTime time_phreeqc;
        time_phreeqc.start();
        chemical_solver_interface_->executeInitialCalculation(
            process_solutions_);
        INFO("[time] Phreeqc took {:g} s.", time_phreeqc.elapsed());
    }

    // All per_process_data_ share the first process.
    bool const is_staggered_coupling =
        !isMonolithicProcess(*per_process_data_[0]);

    if (is_staggered_coupling)
    {
        setCoupledSolutions();
    }

    // Output initial conditions
    {
        const bool output_initial_condition = true;
        outputSolutions(output_initial_condition, 0, start_time_, *output_,
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
    // All per_process_data_ share the first process.
    bool const is_staggered_coupling =
        !isMonolithicProcess(*per_process_data_[0]);

    bool non_equilibrium_initial_residuum_computed = false;
    double t = start_time_;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    NumLib::NonlinearSolverStatus nonlinear_solver_status;

    double dt = computeTimeStepping(0.0, t, accepted_steps, rejected_steps);

    while (t < end_time_)
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
        for (auto& process_data : per_process_data_)
        {
            process_data->process.updateDeactivatedSubdomains(
                t, process_data->process_id);
        }

        if (!non_equilibrium_initial_residuum_computed)
        {
            calculateNonEquilibriumInitialResiduum(
                per_process_data_, process_solutions_, process_solutions_prev_);
            non_equilibrium_initial_residuum_computed = true;
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

        INFO("[time] Time step #{:d} took {:g} s.", timesteps,
             time_timestep.elapsed());

        dt = computeTimeStepping(prev_dt, t, accepted_steps, rejected_steps);

        if (!last_step_rejected_)
        {
            const bool output_initial_condition = false;
            outputSolutions(output_initial_condition, timesteps, t, *output_,
                            &Output::doOutput);
        }

        if (t == end_time_ || t + dt > end_time_ ||
            t + std::numeric_limits<double>::epsilon() > end_time_)
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
                        accepted_steps + rejected_steps, t, *output_,
                        &Output::doOutputLastTimestep);
    }

    return nonlinear_solver_status.error_norms_met;
}

void preTimestepForAllProcesses(
    double const t, double const dt,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<GlobalVector*> const& process_solutions_)
{
    for (auto& process_data : per_process_data)
    {
        auto const process_id = process_data->process_id;
        auto& pcs = process_data->process;
        pcs.preTimestep(process_solutions_, t, dt, process_id);
    }
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

    // All per_process_data_ share the first process.
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
        auto& x = *process_solutions[process_id];
        auto& x_dot = *x_dots[process_id];
        pcs.computeSecondaryVariable(t, dt, x, x_dot, process_id);
        pcs.postTimestep(process_solutions, t, dt, process_id);
    }
}

NumLib::NonlinearSolverStatus TimeLoop::solveUncoupledEquationSystems(
    const double t, const double dt, const std::size_t timestep_id)
{
    preTimestepForAllProcesses(t, dt, per_process_data_, process_solutions_);

    NumLib::NonlinearSolverStatus nonlinear_solver_status;
    for (auto& process_data : per_process_data_)
    {
        auto const process_id = process_data->process_id;
        nonlinear_solver_status = solveMonolithicProcess(
            t, dt, timestep_id, *process_data, process_solutions_,
            process_solutions_prev_, *output_);

        process_data->nonlinear_solver_status = nonlinear_solver_status;
        if (!nonlinear_solver_status.error_norms_met)
        {
            ERR("The nonlinear solver failed in time step #{:d} at t = {:g} s "
                "for process #{:d}.",
                timestep_id, t, process_id);

            if (!process_data->timestepper->canReduceTimestepSize())
            {
                // save unsuccessful solution
                output_->doOutputAlways(process_data->process, process_id,
                                        timestep_id, t, process_solutions_);
                OGS_FATAL(timestepper_cannot_reduce_dt.data());
            }

            return nonlinear_solver_status;
        }
    }

    postTimestepForAllProcesses(t, dt, per_process_data_, process_solutions_,
                                process_solutions_prev_);

    return nonlinear_solver_status;
}

NumLib::NonlinearSolverStatus
TimeLoop::solveCoupledEquationSystemsByStaggeredScheme(
    const double t, const double dt, const std::size_t timestep_id)
{
    // Coupling iteration
    if (global_coupling_max_iterations_ != 0)
    {
        // Set the flag of the first iteration be true.
        for (auto& conv_crit : global_coupling_conv_crit_)
        {
            conv_crit->preFirstIteration();
        }
    }
    auto resetCouplingConvergenceCriteria = [&]() {
        for (auto& conv_crit : global_coupling_conv_crit_)
        {
            conv_crit->reset();
        }
    };

    preTimestepForAllProcesses(t, dt, per_process_data_, process_solutions_);

    NumLib::NonlinearSolverStatus nonlinear_solver_status{false, -1};
    bool coupling_iteration_converged = true;
    for (int global_coupling_iteration = 0;
         global_coupling_iteration < global_coupling_max_iterations_;
         global_coupling_iteration++, resetCouplingConvergenceCriteria())
    {
        // TODO(wenqing): use process name
        coupling_iteration_converged = true;
        int const last_process_id = per_process_data_.size() - 1;
        for (auto& process_data : per_process_data_)
        {
            auto const process_id = process_data->process_id;
            BaseLib::RunTime time_timestep_process;
            time_timestep_process.start();

            CoupledSolutionsForStaggeredScheme coupled_solutions(
                process_solutions_);

            process_data->process.setCoupledSolutionsForStaggeredScheme(
                &coupled_solutions);

            nonlinear_solver_status = solveOneTimeStepOneProcess(
                process_solutions_, process_solutions_prev_, timestep_id, t, dt,
                *process_data, *output_);
            process_data->nonlinear_solver_status = nonlinear_solver_status;

            INFO(
                "[time] Solving process #{:d} took {:g} s in time step #{:d} "
                " coupling iteration #{:d}",
                process_id, time_timestep_process.elapsed(), timestep_id,
                global_coupling_iteration);

            if (!nonlinear_solver_status.error_norms_met)
            {
                ERR("The nonlinear solver failed in time step #{:d} at t = "
                    "{:g} s for process #{:d}.",
                    timestep_id, t, process_id);

                if (!process_data->timestepper->canReduceTimestepSize())
                {
                    // save unsuccessful solution
                    output_->doOutputAlways(process_data->process, process_id,
                                            timestep_id, t, process_solutions_);
                    OGS_FATAL(timestepper_cannot_reduce_dt.data());
                }
                break;
            }

            // Check the convergence of the coupling iteration
            auto& x = *process_solutions_[process_id];
            auto& x_old = *solutions_of_last_cpl_iteration_[process_id];
            if (global_coupling_iteration > 0)
            {
                MathLib::LinAlg::axpy(x_old, -1.0, x);  // save dx to x_old
                if (process_id == last_process_id)
                {
                    INFO(
                        "------- Checking convergence criterion for coupled "
                        "solution  -------");
                    global_coupling_conv_crit_[process_id]->checkDeltaX(x_old,
                                                                        x);
                    coupling_iteration_converged =
                        coupling_iteration_converged &&
                        global_coupling_conv_crit_[process_id]->isSatisfied();
                }
            }
            MathLib::LinAlg::copy(x, x_old);
        }  // end of for (auto& process_data : per_process_data_)

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

    if (chemical_solver_interface_)
    {
        // Sequential non-iterative approach applied here to perform water
        // chemistry calculation followed by resolving component transport
        // process.
        // TODO: move into a global loop to consider both mass balance over
        // space and localized chemical equilibrium between solutes.
        BaseLib::RunTime time_phreeqc;
        time_phreeqc.start();
        chemical_solver_interface_->doWaterChemistryCalculation(
            process_solutions_, dt);
        INFO("[time] Phreeqc took {:g} s.", time_phreeqc.elapsed());
    }

    postTimestepForAllProcesses(t, dt, per_process_data_, process_solutions_,
                                process_solutions_prev_);

    return nonlinear_solver_status;
}

template <typename OutputClass, typename OutputClassMember>
void TimeLoop::outputSolutions(bool const output_initial_condition,
                               unsigned timestep, const double t,
                               OutputClass& output_object,
                               OutputClassMember output_class_member) const
{
    // All per_process_data_ share the first process.
    bool const is_staggered_coupling =
        !isMonolithicProcess(*per_process_data_[0]);

    for (auto& process_data : per_process_data_)
    {
        // If nonlinear solver diverged, the solution has already been
        // saved.
        if (!process_data->nonlinear_solver_status.error_norms_met)
        {
            continue;
        }

        auto const process_id = process_data->process_id;
        auto const& x = *process_solutions_[process_id];
        auto& pcs = process_data->process;

        if (output_initial_condition)
        {
            auto const& ode_sys = *process_data->tdisc_ode_sys;
            // dummy values to handle the time derivative terms more or less
            // correctly, i.e. to ignore them.
            double const t = 0;
            double const dt = 1;
            process_data->time_disc->nextTimestep(t, dt);

            auto& x_dot = NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications(process_id));
            x_dot.setZero();

            pcs.preTimestep(process_solutions_, start_time_, dt, process_id);
            // Update secondary variables, which might be uninitialized, before
            // output.
            pcs.computeSecondaryVariable(start_time_, dt, x, x_dot, process_id);

            NumLib::GlobalVectorProvider::provider.releaseVector(x_dot);
        }
        if (is_staggered_coupling)
        {
            CoupledSolutionsForStaggeredScheme coupled_solutions(
                process_solutions_);

            process_data->process.setCoupledSolutionsForStaggeredScheme(
                &coupled_solutions);
            process_data->process
                .setCoupledTermForTheStaggeredSchemeToLocalAssemblers(
                    process_id);
        }
        (output_object.*output_class_member)(pcs, process_id, timestep, t,
                                             process_solutions_);
    }
}

TimeLoop::~TimeLoop()
{
    for (auto* x : process_solutions_)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }

    for (auto* x : solutions_of_last_cpl_iteration_)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }
}

}  // namespace ProcessLib
