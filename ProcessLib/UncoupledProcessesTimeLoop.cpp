/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UncoupledProcessesTimeLoop.h"

#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/uniqueInsert.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "NumLib/ODESolver/TimeDiscretizationBuilder.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "NumLib/TimeStepping/CreateTimeStepper.h"

#include "MathLib/LinAlg/LinAlg.h"
#include "CoupledSolutionsForStaggeredScheme.h"

std::unique_ptr<ProcessLib::Output> createOutput(
    BaseLib::ConfigTree const& config, std::string const& output_directory)
{
    //! \ogs_file_param{prj__time_loop__output__type}
    config.checkConfigParameter("type", "VTK");
    DBUG("Parse output configuration:");

    return ProcessLib::Output::newInstance(config, output_directory);
}

//! Sets the EquationSystem for the given nonlinear solver,
//! which is Picard or Newton depending on the NLTag.
template <NumLib::NonlinearSolverTag NLTag>
static void setEquationSystem(NumLib::NonlinearSolverBase& nonlinear_solver,
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
static void setEquationSystem(NumLib::NonlinearSolverBase& nonlinear_solver,
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

//! Applies known solutions to the solution vector \c x, transparently
//! for equation systems linearized with either the Picard or Newton method.
template <NumLib::NonlinearSolverTag NLTag>
static void applyKnownSolutions(NumLib::EquationSystem const& eq_sys,
                                GlobalVector& x)
{
    using EqSys = NumLib::NonlinearSystem<NLTag>;
    assert(dynamic_cast<EqSys const*>(&eq_sys) != nullptr);
    auto& eq_sys_ = static_cast<EqSys const&>(eq_sys);

    eq_sys_.applyKnownSolutions(x);
}

//! Applies known solutions to the solution vector \c x, transparently
//! for equation systems linearized with either the Picard or Newton method.
static void applyKnownSolutions(NumLib::EquationSystem const& eq_sys,
                                NumLib::NonlinearSolverTag const nl_tag,
                                GlobalVector& x)
{
    using Tag = NumLib::NonlinearSolverTag;
    switch (nl_tag)
    {
        case Tag::Picard:
            applyKnownSolutions<Tag::Picard>(eq_sys, x);
            break;
        case Tag::Newton:
            applyKnownSolutions<Tag::Newton>(eq_sys, x);
            break;
    }
}

namespace ProcessLib
{
struct SingleProcessData
{
    template <NumLib::NonlinearSolverTag NLTag>
    SingleProcessData(
        std::unique_ptr<NumLib::TimeStepAlgorithm>&& timestepper_,
        NumLib::NonlinearSolver<NLTag>& nonlinear_solver,
        std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit_,
        std::unique_ptr<NumLib::TimeDiscretization>&& time_disc_,
        Process& process_,
        ProcessOutput&& process_output_);

    SingleProcessData(SingleProcessData&& spd);

    std::unique_ptr<NumLib::TimeStepAlgorithm> timestepper;

    //! Flag of skiping time stepping. It is used in the modelling of
    //! coupled processes. If the stepping of any process reaches a steady state
    //! or the ending time, the flag is set to true.
    bool skip_time_stepping = false;

    //! Tag containing the missing type information necessary to cast the
    //! other members of this struct to their concrety types.
    NumLib::NonlinearSolverTag const nonlinear_solver_tag;
    NumLib::NonlinearSolverBase& nonlinear_solver;
    bool nonlinear_solver_converged;
    std::unique_ptr<NumLib::ConvergenceCriterion> conv_crit;

    std::unique_ptr<NumLib::TimeDiscretization> time_disc;
    //! type-erased time-discretized ODE system
    std::unique_ptr<NumLib::EquationSystem> tdisc_ode_sys;
    //! cast of \c tdisc_ode_sys to NumLib::InternalMatrixStorage
    NumLib::InternalMatrixStorage* mat_strg = nullptr;

    Process& process;
    ProcessOutput process_output;
};

template <NumLib::NonlinearSolverTag NLTag>
SingleProcessData::SingleProcessData(
    std::unique_ptr<NumLib::TimeStepAlgorithm>&& timestepper_,
    NumLib::NonlinearSolver<NLTag>& nonlinear_solver,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit_,
    std::unique_ptr<NumLib::TimeDiscretization>&& time_disc_,
    Process& process_,
    ProcessOutput&& process_output_)
    : timestepper(std::move(timestepper_)),
      nonlinear_solver_tag(NLTag),
      nonlinear_solver(nonlinear_solver),
      nonlinear_solver_converged(true),
      conv_crit(std::move(conv_crit_)),
      time_disc(std::move(time_disc_)),
      process(process_),
      process_output(std::move(process_output_))
{
}

SingleProcessData::SingleProcessData(SingleProcessData&& spd)
    : timestepper(std::move(spd.timestepper)),
      nonlinear_solver_tag(spd.nonlinear_solver_tag),
      nonlinear_solver(spd.nonlinear_solver),
      nonlinear_solver_converged(spd.nonlinear_solver_converged),
      conv_crit(std::move(spd.conv_crit)),
      time_disc(std::move(spd.time_disc)),
      tdisc_ode_sys(std::move(spd.tdisc_ode_sys)),
      mat_strg(spd.mat_strg),
      process(spd.process),
      process_output(std::move(spd.process_output))
{
    spd.mat_strg = nullptr;
}

template <NumLib::ODESystemTag ODETag>
void setTimeDiscretizedODESystem(
    SingleProcessData& spd,
    NumLib::ODESystem<ODETag, NumLib::NonlinearSolverTag::Picard>& ode_sys)
{
    using Tag = NumLib::NonlinearSolverTag;
    // A concrete Picard solver
    using NonlinearSolverPicard = NumLib::NonlinearSolver<Tag::Picard>;
    // A concrete Newton solver
    using NonlinearSolverNewton = NumLib::NonlinearSolver<Tag::Newton>;

    if (dynamic_cast<NonlinearSolverPicard*>(&spd.nonlinear_solver))
    {
        // The Picard solver can also work with a Newton-ready ODE,
        // because the Newton ODESystem derives from the Picard ODESystem.
        // So no further checks are needed here.

        spd.tdisc_ode_sys = std::make_unique<
            NumLib::TimeDiscretizedODESystem<ODETag, Tag::Picard>>(
            ode_sys, *spd.time_disc);
    }
    else if (dynamic_cast<NonlinearSolverNewton*>(&spd.nonlinear_solver))
    {
        // The Newton-Raphson method needs a Newton-ready ODE.

        using ODENewton = NumLib::ODESystem<ODETag, Tag::Newton>;
        if (auto* ode_newton = dynamic_cast<ODENewton*>(&ode_sys))
        {
            spd.tdisc_ode_sys = std::make_unique<
                NumLib::TimeDiscretizedODESystem<ODETag, Tag::Newton>>(
                *ode_newton, *spd.time_disc);
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

    spd.mat_strg =
        dynamic_cast<NumLib::InternalMatrixStorage*>(spd.tdisc_ode_sys.get());
}

void setTimeDiscretizedODESystem(SingleProcessData& spd)
{
    setTimeDiscretizedODESystem(spd, spd.process);
}

std::unique_ptr<SingleProcessData> makeSingleProcessData(
    std::unique_ptr<NumLib::TimeStepAlgorithm>&& timestepper,
    NumLib::NonlinearSolverBase& nonlinear_solver,
    Process& process,
    std::unique_ptr<NumLib::TimeDiscretization>&& time_disc,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit,
    ProcessOutput&& process_output)
{
    using Tag = NumLib::NonlinearSolverTag;

    if (auto* nonlinear_solver_picard =
            dynamic_cast<NumLib::NonlinearSolver<Tag::Picard>*>(
                &nonlinear_solver))
    {
        return std::make_unique<SingleProcessData>(
            std::move(timestepper), *nonlinear_solver_picard,
            std::move(conv_crit), std::move(time_disc), process,
            std::move(process_output));
    }
    if (auto* nonlinear_solver_newton =
            dynamic_cast<NumLib::NonlinearSolver<Tag::Newton>*>(
                &nonlinear_solver))
    {
        return std::make_unique<SingleProcessData>(
            std::move(timestepper), *nonlinear_solver_newton,
            std::move(conv_crit), std::move(time_disc), process,
            std::move(process_output));
    }

    OGS_FATAL("Encountered unknown nonlinear solver type. Aborting");
}

std::vector<std::unique_ptr<SingleProcessData>> createPerProcessData(
    BaseLib::ConfigTree const& config,
    const std::map<std::string, std::unique_ptr<Process>>& processes,
    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>> const&
        nonlinear_solvers)
{
    std::vector<std::unique_ptr<SingleProcessData>> per_process_data;

    //! \ogs_file_param{prj__time_loop__processes__process}
    for (auto pcs_config : config.getConfigSubtreeList("process"))
    {
        //! \ogs_file_attr{prj__time_loop__processes__process__ref}
        auto const pcs_name = pcs_config.getConfigAttribute<std::string>("ref");
        auto& pcs = *BaseLib::getOrError(
            processes, pcs_name,
            "A process with the given name has not been defined.");

        auto const nl_slv_name =
            //! \ogs_file_param{prj__time_loop__processes__process__nonlinear_solver}
            pcs_config.getConfigParameter<std::string>("nonlinear_solver");
        auto& nl_slv = *BaseLib::getOrError(
            nonlinear_solvers, nl_slv_name,
            "A nonlinear solver with the given name has not been defined.");

        auto time_disc = NumLib::createTimeDiscretization(
            //! \ogs_file_param{prj__time_loop__processes__process__time_discretization}
            pcs_config.getConfigSubtree("time_discretization"));

        auto timestepper = NumLib::createTimeStepper(
            //! \ogs_file_param{prj__time_loop__processes__process__time_stepping}
            pcs_config.getConfigSubtree("time_stepping"));

        auto conv_crit = NumLib::createConvergenceCriterion(
            //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion}
            pcs_config.getConfigSubtree("convergence_criterion"));

        //! \ogs_file_param{prj__time_loop__processes__process__output}
        ProcessOutput process_output{pcs_config.getConfigSubtree("output")};

        per_process_data.emplace_back(makeSingleProcessData(
            std::move(timestepper), nl_slv, pcs, std::move(time_disc),
            std::move(conv_crit), std::move(process_output)));
    }

    if (per_process_data.size() != processes.size())
    {
        if (processes.size() > 1)
        {
            OGS_FATAL(
                "Some processes have not been configured to be solved by this "
                " time loop.");
        }
        else
        {
            INFO(
                "The equations of the coupled processes will be solved by the "
                "staggered scheme.")
        }
    }

    return per_process_data;
}

std::unique_ptr<UncoupledProcessesTimeLoop> createUncoupledProcessesTimeLoop(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    const std::map<std::string, std::unique_ptr<Process>>& processes,
    const std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>>&
        nonlinear_solvers)
{
    auto const& coupling_config
        //! \ogs_file_param{prj__time_loop__global_process_coupling}
        = config.getConfigSubtreeOptional("global_process_coupling");

    std::unique_ptr<NumLib::ConvergenceCriterion> coupling_conv_crit = nullptr;
    unsigned max_coupling_iterations = 1;
    if (coupling_config)
    {
        max_coupling_iterations
            //! \ogs_file_param{prj__time_loop__global_process_coupling__max_iter}
            = coupling_config->getConfigParameter<unsigned>("max_iter");
        coupling_conv_crit = NumLib::createConvergenceCriterion(
            //! \ogs_file_param{prj__time_loop__global_process_coupling__convergence_criterion}
            coupling_config->getConfigSubtree("convergence_criterion"));
    }

    auto output =
        //! \ogs_file_param{prj__time_loop__output}
        createOutput(config.getConfigSubtree("output"), output_directory);

    auto per_process_data = createPerProcessData(
        //! \ogs_file_param{prj__time_loop__processes}
        config.getConfigSubtree("processes"), processes, nonlinear_solvers);

    const auto minmax_iter = std::minmax_element(
        per_process_data.begin(),
        per_process_data.end(),
        [](std::unique_ptr<SingleProcessData> const& a,
           std::unique_ptr<SingleProcessData> const& b) {
            return (a->timestepper->end() < b->timestepper->end());
        });
    const double start_time =
        per_process_data[minmax_iter.first - per_process_data.begin()]
            ->timestepper->begin();
    const double end_time =
        per_process_data[minmax_iter.second - per_process_data.begin()]
            ->timestepper->end();

    return std::make_unique<UncoupledProcessesTimeLoop>(
        std::move(output), std::move(per_process_data), max_coupling_iterations,
        std::move(coupling_conv_crit), start_time, end_time);
}

std::vector<GlobalVector*> setInitialConditions(
    double const t0,
    std::vector<std::unique_ptr<SingleProcessData>> const& per_process_data)
{
    std::vector<GlobalVector*> process_solutions;

    unsigned pcs_idx = 0;
    for (auto& spd : per_process_data)
    {
        auto& pcs = spd->process;
        auto& time_disc = *spd->time_disc;

        auto& ode_sys = *spd->tdisc_ode_sys;
        auto const nl_tag = spd->nonlinear_solver_tag;

        // append a solution vector of suitable size
        process_solutions.emplace_back(
            &NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications()));

        auto& x0 = *process_solutions[pcs_idx];
        pcs.setInitialConditions(pcs_idx, t0, x0);
        MathLib::LinAlg::finalizeAssembly(x0);

        time_disc.setInitialState(t0, x0);  // push IC

        if (time_disc.needsPreload())
        {
            auto& nonlinear_solver = spd->nonlinear_solver;
            auto& mat_strg = *spd->mat_strg;
            auto& conv_crit = *spd->conv_crit;

            setEquationSystem(nonlinear_solver, ode_sys, conv_crit, nl_tag);
            nonlinear_solver.assemble(x0);
            time_disc.pushState(
                t0, x0,
                mat_strg);  // TODO: that might do duplicate work
        }

        ++pcs_idx;
    }

    return process_solutions;
}

bool solveOneTimeStepOneProcess(const unsigned process_index,
                                GlobalVector& x, std::size_t const timestep,
                                double const t, double const delta_t,
                                SingleProcessData const& process_data,
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

    applyKnownSolutions(ode_sys, nl_tag, x);

    auto const post_iteration_callback = [&](unsigned iteration,
                                             GlobalVector const& x) {
        output_control.doOutputNonlinearIteration(
            process, process_index, process_data.process_output, timestep, t,
            x, iteration);
    };

    bool nonlinear_solver_succeeded =
        nonlinear_solver.solve(x, post_iteration_callback);

    return nonlinear_solver_succeeded;
}

UncoupledProcessesTimeLoop::UncoupledProcessesTimeLoop(
    std::unique_ptr<Output>&& output,
    std::vector<std::unique_ptr<SingleProcessData>>&& per_process_data,
    const unsigned global_coupling_max_iterations,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& global_coupling_conv_crit,
    const double start_time, const double end_time)
    : _output(std::move(output)),
      _per_process_data(std::move(per_process_data)),
      _start_time(start_time),
      _end_time(end_time),
      _global_coupling_max_iterations(global_coupling_max_iterations),
      _global_coupling_conv_crit(std::move(global_coupling_conv_crit))
{
}

bool UncoupledProcessesTimeLoop::setCoupledSolutions()
{
    // All _per_process_data share one process
    const bool use_monolithic_scheme =
        _per_process_data[0]->process.isMonolithicSchemeUsed();
    if (use_monolithic_scheme)
    {
        return false;
    }

    _solutions_of_coupled_processes.reserve(_per_process_data.size());
    for (unsigned pcs_idx = 0; pcs_idx < _per_process_data.size(); pcs_idx++)
    {
        auto const& x = *_process_solutions[pcs_idx];
        _solutions_of_coupled_processes.emplace_back(x);

        // Create a vector to store the solution of the last coupling iteration
        auto& x0 = NumLib::GlobalVectorProvider::provider.getVector(x);
        MathLib::LinAlg::copy(x, x0);

        // append a solution vector of suitable size
        _solutions_of_last_cpl_iteration.emplace_back(&x0);
    }

    return true;  // use staggered scheme.
}

double UncoupledProcessesTimeLoop::computeTimeStepping(
    const double prev_dt, double& t, std::size_t& accepted_steps,
    std::size_t& rejected_steps)
{
    bool all_process_steps_accepted = true;
    // Get minimum time step size among step sizes of all processes.
    double dt = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < _per_process_data.size(); i++)
    {
        auto& ppd = *_per_process_data[i];
        const auto& timestepper = ppd.timestepper;
        if (t > timestepper->end())
        {
            // skip the process that already reaches the ending time.
            ppd.skip_time_stepping = true;
            continue;
        }

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

        if (!ppd.nonlinear_solver_converged)
        {
            timestepper->setAcceptedOrNot(false);
        }

        if (!timestepper->next(solution_error) &&
            // In case of FixedTimeStepping, which makes timestepper->next(...)
            // return false when the ending time is reached.
            t + std::numeric_limits<double>::epsilon() < timestepper->end())
        {
            // Not all processes have accepted steps.
            all_process_steps_accepted = false;
        }

        if (!ppd.nonlinear_solver_converged)
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
        auto& mat_strg = *ppd.mat_strg;
        auto& x = *_process_solutions[i];
        if (all_process_steps_accepted)
        {
            time_disc->pushState(t, x, mat_strg);
        }
        else
        {
            if (t < _end_time)
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
            if (t < _end_time)
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
bool UncoupledProcessesTimeLoop::loop()
{
    // initialize output, convergence criterion, etc.
    {
        unsigned pcs_idx = 0;
        for (auto& spd : _per_process_data)
        {
            auto& pcs = spd->process;
            _output->addProcess(pcs, pcs_idx);

            setTimeDiscretizedODESystem(*spd);

            if (auto* conv_crit =
                    dynamic_cast<NumLib::ConvergenceCriterionPerComponent*>(
                        spd->conv_crit.get()))
            {
                conv_crit->setDOFTable(pcs.getDOFTable(), pcs.getMesh());
            }

            ++pcs_idx;
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
    bool nonlinear_solver_succeeded = true;

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

        if (is_staggered_coupling)
        {
            nonlinear_solver_succeeded =
                solveCoupledEquationSystemsByStaggeredScheme(t, dt, timesteps);
        }
        else
        {
            nonlinear_solver_succeeded =
                solveUncoupledEquationSystems(t, dt, timesteps);
        }

        INFO("[time] Time step #%u took %g s.", timesteps,
             time_timestep.elapsed());

        dt = computeTimeStepping(prev_dt, t, accepted_steps, rejected_steps);

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

        // If this step was rejected twice with the same time step size, jump
        // out this function directly with a failure flag and let the main
        // function terminate the program.
        if (std::abs(dt - prev_dt) < std::numeric_limits<double>::min() &&
            _last_step_rejected)
        {
            ALERT(
                "\tTime step %u is rejected and the new computed"
                " step size is the same as\n"
                "\tthat was just used.\n"
                "\tSuggest to adjust the parameters of the time"
                " stepper or try other time stepper.\n"
                "\tThe program will stop.",
                timesteps);
            // save unsuccessful solution
            const bool output_initial_condition = false;
            outputSolutions(output_initial_condition, is_staggered_coupling,
                            timesteps, t, *_output, &Output::doOutputAlways);
            return false;
        }
    }

    INFO(
        "The whole computation of the time stepping took %u steps, in which\n"
        "\t the accepted steps are %u, and the rejected steps are %u.\n",
        accepted_steps + rejected_steps, accepted_steps, rejected_steps)

    // output last time step
    if (nonlinear_solver_succeeded)
    {
        const bool output_initial_condition = false;
        outputSolutions(output_initial_condition, is_staggered_coupling,
                        accepted_steps + rejected_steps, t, *_output,
                        &Output::doOutputLastTimestep);
    }

    return nonlinear_solver_succeeded;
}

static std::string nonlinear_fixed_dt_fails_info =
    "Nonlinear solver fails. Because of the time stepper"
    " of FixedTimeStepping is used, the program has to be"
    " terminated ";

bool UncoupledProcessesTimeLoop::solveUncoupledEquationSystems(
    const double t, const double dt, const std::size_t timestep_id)
{
    // TODO(wenqing): use process name
    unsigned pcs_idx = 0;
    for (auto& spd : _per_process_data)
    {
        if (spd->skip_time_stepping)
        {
            INFO("Process %u is skipped in the time stepping.", pcs_idx);
            ++pcs_idx;
            continue;
        }

        BaseLib::RunTime time_timestep_process;
        time_timestep_process.start();

        auto& x = *_process_solutions[pcs_idx];
        auto& pcs = spd->process;
        pcs.preTimestep(x, t, dt, pcs_idx);

        const auto nonlinear_solver_succeeded =
            solveOneTimeStepOneProcess(pcs_idx, x, timestep_id, t, dt, *spd, *_output);
        spd->nonlinear_solver_converged = nonlinear_solver_succeeded;
        pcs.postTimestep(x);
        pcs.computeSecondaryVariable(t, x);

        INFO("[time] Solving process #%u took %g s in time step #%u ", pcs_idx,
             time_timestep_process.elapsed(), timestep_id);

        if (!nonlinear_solver_succeeded)
        {
            ERR("The nonlinear solver failed in time step #%u at t = %g "
                "s for process #%u.",
                timestep_id, t, pcs_idx);

            // save unsuccessful solution
            _output->doOutputAlways(pcs, pcs_idx, spd->process_output, timestep_id, t,
                                    x);

            if (!spd->timestepper->isSolutionErrorComputationNeeded())
            {
                OGS_FATAL(nonlinear_fixed_dt_fails_info.data());
            }

            return false;
        }

        _output->doOutput(pcs, pcs_idx, spd->process_output, timestep_id, t, x);

        ++pcs_idx;
    }  // end of for (auto& spd : _per_process_data)

    return true;
}

bool UncoupledProcessesTimeLoop::solveCoupledEquationSystemsByStaggeredScheme(
    const double t, const double dt, const std::size_t timestep_id)
{
    // Coupling iteration
    if (_global_coupling_max_iterations != 0)
    {
        // Set the flag of the first iteration be true.
        _global_coupling_conv_crit->preFirstIteration();
    }
    bool coupling_iteration_converged = true;
    for (unsigned global_coupling_iteration = 0;
         global_coupling_iteration < _global_coupling_max_iterations;
         global_coupling_iteration++, _global_coupling_conv_crit->reset())
    {
        // TODO(wenqing): use process name
        bool nonlinear_solver_succeeded = true;
        coupling_iteration_converged = true;
        unsigned pcs_idx = 0;
        for (auto& spd : _per_process_data)
        {
            if (spd->skip_time_stepping)
            {
                INFO("Process %u is skipped in the time stepping.", pcs_idx);
                ++pcs_idx;
                continue;
            }

            BaseLib::RunTime time_timestep_process;
            time_timestep_process.start();

            auto& x = *_process_solutions[pcs_idx];
            if (global_coupling_iteration == 0)
            {
                // Copy the solution of the previous time step to a vector that
                // belongs to process. For some problems, both of the current
                // solution and the solution of the previous time step are
                // required for the coupling computation.
                spd->process.preTimestep(x, t, dt, pcs_idx);
            }

            CoupledSolutionsForStaggeredScheme coupled_solutions(
                _solutions_of_coupled_processes, dt, pcs_idx);

            spd->process.setCoupledSolutionsForStaggeredScheme(
                &coupled_solutions);

            const auto nonlinear_solver_succeeded = solveOneTimeStepOneProcess(
                pcs_idx, x, timestep_id, t, dt, *spd, *_output);
            spd->nonlinear_solver_converged = nonlinear_solver_succeeded;

            INFO(
                "[time] Solving process #%u took %g s in time step #%u "
                " coupling iteration #%u",
                pcs_idx, time_timestep_process.elapsed(), timestep_id,
                global_coupling_iteration);

            if (!nonlinear_solver_succeeded)
            {
                ERR("The nonlinear solver failed in time step #%u at t = %g "
                    "s"
                    " for process #%u.",
                    timestep_id, t, pcs_idx);

                // save unsuccessful solution
                _output->doOutputAlways(spd->process, pcs_idx,
                                        spd->process_output,
                                        timestep_id, t, x);

                if (!spd->timestepper->isSolutionErrorComputationNeeded())
                {
                    OGS_FATAL(nonlinear_fixed_dt_fails_info.data());
                }

                break;
            }

            // Check the convergence of the coupling iteration
            auto& x_old = *_solutions_of_last_cpl_iteration[pcs_idx];
            MathLib::LinAlg::axpy(x_old, -1.0, x);
            _global_coupling_conv_crit->checkResidual(x_old);
            coupling_iteration_converged =
                coupling_iteration_converged &&
                _global_coupling_conv_crit->isSatisfied();

            MathLib::LinAlg::copy(x, x_old);

            if (coupling_iteration_converged)
            {
                break;
            }

            ++pcs_idx;
        }  // end of for (auto& spd : _per_process_data)

        if (coupling_iteration_converged)
        {
            break;
        }

        if (!nonlinear_solver_succeeded)
        {
            return false;
        }
    }

    if (!coupling_iteration_converged)
    {
        WARN(
            "The coupling iterations reaches its maximum number in time step"
            "#%u at t = %g s",
            timestep_id, t);
    }

    unsigned pcs_idx = 0;
    for (auto& spd : _per_process_data)
    {
        if (spd->skip_time_stepping)
        {
            ++pcs_idx;
            continue;
        }
        auto& pcs = spd->process;
        auto& x = *_process_solutions[pcs_idx];
        pcs.postTimestep(x);
        pcs.computeSecondaryVariable(t, x);

        ++pcs_idx;
    }

    {
        const bool output_initial_condition = false;
        const bool is_staggered_coupling = true;
        outputSolutions(output_initial_condition, is_staggered_coupling,
                        timestep_id, t, *_output, &Output::doOutput);
    }

    return true;
}

template <typename OutputClass, typename OutputClassMember>
void UncoupledProcessesTimeLoop::outputSolutions(
    bool const output_initial_condition, bool const is_staggered_coupling,
    unsigned timestep, const double t, OutputClass& output_object,
    OutputClassMember output_class_member) const
{
    unsigned pcs_idx = 0;
    for (auto& spd : _per_process_data)
    {
        auto& pcs = spd->process;
        // If nonlinear solver diverged, the solution has already been
        // saved.
        if ((!spd->nonlinear_solver_converged) || spd->skip_time_stepping)
        {
            ++pcs_idx;
            continue;
        }

        auto const& x = *_process_solutions[pcs_idx];

        if (output_initial_condition)
        {
            pcs.preTimestep(x, _start_time,
                            spd->timestepper->getTimeStep().dt(), pcs_idx);
        }
        if (is_staggered_coupling)
        {
            CoupledSolutionsForStaggeredScheme coupled_solutions(
                _solutions_of_coupled_processes, 0.0, pcs_idx);

            spd->process.setCoupledSolutionsForStaggeredScheme(
                &coupled_solutions);
            spd->process.setCoupledTermForTheStaggeredSchemeToLocalAssemblers();
            (output_object.*output_class_member)(pcs, pcs_idx,
                                                 spd->process_output,
                                                 timestep, t, x);
        }
        else
        {
            (output_object.*output_class_member)(
                pcs, pcs_idx, spd->process_output, timestep, t, x);
        }

        ++pcs_idx;
    }
}

UncoupledProcessesTimeLoop::~UncoupledProcessesTimeLoop()
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
