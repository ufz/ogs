/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UncoupledProcessesTimeLoop.h"
#include "BaseLib/uniqueInsert.h"
#include "BaseLib/RunTime.h"
#include "NumLib/ODESolver/TimeDiscretizationBuilder.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

#include "MathLib/LinAlg/LinAlg.h"

std::unique_ptr<NumLib::ITimeStepAlgorithm> createTimeStepper(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__time_stepping__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    std::unique_ptr<NumLib::ITimeStepAlgorithm> timestepper;

    if (type == "SingleStep")
    {
        //! \ogs_file_param_special{prj__time_loop__time_stepping__SingleStep}
        config.ignoreConfigParameter("type");
        timestepper.reset(new NumLib::FixedTimeStepping(0.0, 1.0, 1.0));
    }
    else if (type == "FixedTimeStepping")
    {
        timestepper = NumLib::FixedTimeStepping::newInstance(config);
    }
    else
    {
        OGS_FATAL("Unknown timestepper type: `%s'.", type.c_str());
    }

    return timestepper;
}

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
        NumLib::NonlinearSolver<NLTag>& nonlinear_solver,
        std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit_,
        std::unique_ptr<NumLib::TimeDiscretization>&& time_disc_,
        Process& process_,
        std::map<ProcessType, Process const&>&& coupled_processes_,
        ProcessOutput&& process_output_);

    SingleProcessData(SingleProcessData&& spd);

    //! Tag containing the missing type information necessary to cast the
    //! other members of this struct to their concrety types.
    NumLib::NonlinearSolverTag const nonlinear_solver_tag;
    NumLib::NonlinearSolverBase& nonlinear_solver;
    std::unique_ptr<NumLib::ConvergenceCriterion> conv_crit;

    std::unique_ptr<NumLib::TimeDiscretization> time_disc;
    //! type-erased time-discretized ODE system
    std::unique_ptr<NumLib::EquationSystem> tdisc_ode_sys;
    //! cast of \c tdisc_ode_sys to NumLib::InternalMatrixStorage
    NumLib::InternalMatrixStorage* mat_strg = nullptr;

    Process& process;
    /// Coupled processes.
    std::map<ProcessType, Process const&> const coupled_processes;
    ProcessOutput process_output;
};

template <NumLib::NonlinearSolverTag NLTag>
SingleProcessData::SingleProcessData(
    NumLib::NonlinearSolver<NLTag>& nonlinear_solver,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit_,
    std::unique_ptr<NumLib::TimeDiscretization>&& time_disc_,
    Process& process_,
    std::map<ProcessType, Process const&>&& coupled_processes_,
    ProcessOutput&& process_output_)
    : nonlinear_solver_tag(NLTag),
      nonlinear_solver(nonlinear_solver),
      conv_crit(std::move(conv_crit_)),
      time_disc(std::move(time_disc_)),
      process(process_),
      coupled_processes(coupled_processes_),
      process_output(std::move(process_output_))
{
}

SingleProcessData::SingleProcessData(SingleProcessData&& spd)
    : nonlinear_solver_tag(spd.nonlinear_solver_tag),
      nonlinear_solver(spd.nonlinear_solver),
      conv_crit(std::move(spd.conv_crit)),
      time_disc(std::move(spd.time_disc)),
      tdisc_ode_sys(std::move(spd.tdisc_ode_sys)),
      mat_strg(spd.mat_strg),
      process(spd.process),
      coupled_processes(spd.coupled_processes),
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

        spd.tdisc_ode_sys.reset(
            new NumLib::TimeDiscretizedODESystem<ODETag, Tag::Picard>(
                ode_sys, *spd.time_disc));
    }
    else if (dynamic_cast<NonlinearSolverNewton*>(&spd.nonlinear_solver))
    {
        // The Newton-Raphson method needs a Newton-ready ODE.

        using ODENewton = NumLib::ODESystem<ODETag, Tag::Newton>;
        if (auto* ode_newton = dynamic_cast<ODENewton*>(&ode_sys))
        {
            spd.tdisc_ode_sys.reset(
                new NumLib::TimeDiscretizedODESystem<ODETag, Tag::Newton>(
                    *ode_newton, *spd.time_disc));
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
    NumLib::NonlinearSolverBase& nonlinear_solver,
    Process& process,
    std::unique_ptr<NumLib::TimeDiscretization>&& time_disc,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit,
    std::map<ProcessType, Process const&>&& coupled_processes,
    ProcessOutput&& process_output)
{
    using Tag = NumLib::NonlinearSolverTag;

    if (auto* nonlinear_solver_picard =
            dynamic_cast<NumLib::NonlinearSolver<Tag::Picard>*>(
                &nonlinear_solver))
    {
        return std::unique_ptr<SingleProcessData>{new SingleProcessData{
            *nonlinear_solver_picard, std::move(conv_crit),
            std::move(time_disc), process, std::move(coupled_processes),
            std::move(process_output)}};
    }
    else if (auto* nonlinear_solver_newton =
                 dynamic_cast<NumLib::NonlinearSolver<Tag::Newton>*>(
                     &nonlinear_solver))
    {
        return std::unique_ptr<SingleProcessData>{new SingleProcessData{
            *nonlinear_solver_newton, std::move(conv_crit),
            std::move(time_disc), process, std::move(coupled_processes),
            std::move(process_output)}};
    }
    else
    {
        OGS_FATAL("Encountered unknown nonlinear solver type. Aborting");
    }
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

        auto conv_crit = NumLib::createConvergenceCriterion(
            //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion}
            pcs_config.getConfigSubtree("convergence_criterion"));

        auto const& cpl_pcses
            //! \ogs_file_param{prj__time_loop__processes__process__coupled_processes}
            = pcs_config.getConfigSubtreeOptional("coupled_processes");
        std::map<ProcessType, Process const&> coupled_processes;
        if (cpl_pcses)
        {
            for (
                auto const cpl_pcs_name :
                //! \ogs_file_param{prj__time_loop__processes__process__coupled_processes__coupled_process}
                cpl_pcses->getConfigParameterList<std::string>(
                    "coupled_process"))
            {
                auto& cpl_pcs_ptr = BaseLib::getOrError(
                    processes, cpl_pcs_name,
                    "A process with the given name has not been defined.");

                auto const cpl_pcs = cpl_pcs_ptr.get();
                auto const inserted = coupled_processes.emplace(
                    cpl_pcs_ptr->getProcessType(), *cpl_pcs);
                if (!inserted.second)
                {  // insertion failed, i.e., key already exists
                    OGS_FATAL("Coupled process `%s' already exists.",
                              cpl_pcs_name.data());
                }
            }
        }

        //! \ogs_file_param{prj__time_loop__processes__process__output}
        ProcessOutput process_output{pcs_config.getConfigSubtree("output")};

        per_process_data.emplace_back(makeSingleProcessData(
            nl_slv, pcs, std::move(time_disc), std::move(conv_crit),
            std::move(coupled_processes), std::move(process_output)));
    }

    if (per_process_data.size() != processes.size())
        OGS_FATAL(
            "Some processes have not been configured to be solved by this time "
            "time loop.");

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
            //! \ogs_file_param{prj__time_loop__global_process_coupling__max_iteration}
            = coupling_config->getConfigParameter<unsigned>("max_iteration");
        coupling_conv_crit = NumLib::createConvergenceCriterion(
            //! \ogs_file_param{prj__time_loop__global_process_coupling__max_iteration__convergence_criterion}
            coupling_config->getConfigSubtree("convergence_criterion"));
    }

    auto timestepper =
        //! \ogs_file_param{prj__time_loop__time_stepping}
        createTimeStepper(config.getConfigSubtree("time_stepping"));

    auto output =
        //! \ogs_file_param{prj__time_loop__output}
        createOutput(config.getConfigSubtree("output"), output_directory);

    auto per_process_data = createPerProcessData(
        //! \ogs_file_param{prj__time_loop__processes}
        config.getConfigSubtree("processes"), processes, nonlinear_solvers);

    return std::unique_ptr<UncoupledProcessesTimeLoop>{
        new UncoupledProcessesTimeLoop{
            std::move(timestepper), std::move(output),
            std::move(per_process_data), max_coupling_iterations,
            std::move(coupling_conv_crit)}};
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
        pcs.setInitialConditions(t0, x0);
        MathLib::LinAlg::finalizeAssembly(x0);

        time_disc.setInitialState(t0, x0);  // push IC

        if (time_disc.needsPreload())
        {
            auto& nonlinear_solver = spd->nonlinear_solver;
            auto& mat_strg = *spd->mat_strg;
            auto& conv_crit = *spd->conv_crit;

            setEquationSystem(nonlinear_solver, ode_sys, conv_crit, nl_tag);
            nonlinear_solver.assemble(
                x0, ProcessLib::createVoidStaggeredCouplingTerm());
            time_disc.pushState(
                t0, x0, mat_strg);  // TODO: that might do duplicate work
        }

        ++pcs_idx;
    }

    return process_solutions;
}

bool solveOneTimeStepOneProcess(GlobalVector& x, std::size_t const timestep,
                                double const t, double const delta_t,
                                SingleProcessData& process_data,
                                StaggeredCouplingTerm const& coupled_term,
                                Output const& output_control)
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
            process, process_data.process_output, timestep, t, x, iteration);
    };

    bool nonlinear_solver_succeeded =
        nonlinear_solver.solve(x, coupled_term, post_iteration_callback);

    auto& mat_strg = *process_data.mat_strg;
    time_disc.pushState(t, x, mat_strg);

    return nonlinear_solver_succeeded;
}

UncoupledProcessesTimeLoop::UncoupledProcessesTimeLoop(
    std::unique_ptr<NumLib::ITimeStepAlgorithm>&& timestepper,
    std::unique_ptr<Output>&& output,
    std::vector<std::unique_ptr<SingleProcessData>>&& per_process_data,
    const unsigned global_coupling_max_iterations,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& glb_coupling_conv_crit)
    : _timestepper{std::move(timestepper)},
      _output(std::move(output)),
      _per_process_data(std::move(per_process_data)),
      _global_coupling_max_iterations(global_coupling_max_iterations),
      _global_coupling_conv_crit(std::move(glb_coupling_conv_crit))
{
}

/**
 *  This function fills the vector of solutions of coupled processes of
 *  processes, _solutions_of_coupled_processes, and initializes the vector of
 *  solutions of the previous coupling iteration,
 *  _solutions_of_last_cpl_iteration.
 *
 *  \return a boolean value as a flag to indicate there should be a coupling
 *          among process or not.
 */
bool UncoupledProcessesTimeLoop::setCoupledSolutions()
{
    if (!_global_coupling_conv_crit && _global_coupling_max_iterations == 1)
        return false;

    unsigned pcs_idx = 0;
    _solutions_of_coupled_processes.reserve(_per_process_data.size());
    for (auto& spd : _per_process_data)
    {
        BaseLib::RunTime time_timestep_process;
        time_timestep_process.start();

        auto const& coupled_processes = spd->coupled_processes;
        std::map<ProcessType, GlobalVector const&> coupled_xs;
        auto it = coupled_processes.begin();
        while (it != coupled_processes.end())
        {
            ProcessLib::Process const& coupled_process = it->second;
            auto const found_item = std::find_if(
                _per_process_data.begin(),
                _per_process_data.end(),
                [&coupled_process](
                    std::unique_ptr<SingleProcessData> const& item) {
                    return coupled_process.getProcessType() ==
                           item->process.getProcessType();
                });

            if (found_item != _per_process_data.end())
            {
                // Id of the coupled process:
                const std::size_t c_id = found_item - _per_process_data.begin();

                BaseLib::insertIfKeyUniqueElseError(
                    coupled_xs, it->first, *_process_solutions[c_id],
                    "global_coupled_x");
            }
            it++;
        }
        _solutions_of_coupled_processes.emplace_back(coupled_xs);

        const auto x = _process_solutions[pcs_idx];

        auto x1 = &NumLib::GlobalVectorProvider::provider.getVector(*x);
        MathLib::LinAlg::copy(*x, *x1);

        // append a solution vector of suitable size
        _solutions_of_last_cpl_iteration.emplace_back(x1);

        ++pcs_idx;
    }

    return true;
}

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

    auto const t0 = _timestepper->getTimeStep().current();  // time of the IC

    // init solution storage
    _process_solutions = setInitialConditions(t0, _per_process_data);

    // output initial conditions
    {
        unsigned pcs_idx = 0;
        for (auto& spd : _per_process_data)
        {
            auto& pcs = spd->process;
            auto const& x0 = *_process_solutions[pcs_idx];

            pcs.preTimestep(x0, t0, _timestepper->getTimeStep().dt());
            _output->doOutput(pcs, spd->process_output, 0, t0, x0);
            ++pcs_idx;
        }
    }

    const bool is_staggered_coupling = setCoupledSolutions();

    double t = t0;
    std::size_t timestep = 1;  // the first timestep really is number one
    bool nonlinear_solver_succeeded = true;

    while (_timestepper->next())
    {
        BaseLib::RunTime time_timestep;
        time_timestep.start();

        auto const ts = _timestepper->getTimeStep();
        auto const delta_t = ts.dt();
        t = ts.current();
        timestep = ts.steps();

        INFO("=== timestep #%u (t=%gs, dt=%gs) ==============================",
             timestep, t, delta_t);

        if (is_staggered_coupling)
            nonlinear_solver_succeeded =
                solveCoupledEquationSystemsByStaggeredScheme(t, delta_t,
                                                             timestep);
        else
            nonlinear_solver_succeeded =
                solveUncoupledEquationSystems(t, delta_t, timestep);

        INFO("[time] Timestep #%u took %g s.", timestep,
             time_timestep.elapsed());

        if (!nonlinear_solver_succeeded)
            break;
    }

    // output last timestep
    if (nonlinear_solver_succeeded)
    {
        unsigned pcs_idx = 0;
        for (auto& spd : _per_process_data)
        {
            auto& pcs = spd->process;
            auto const& x = *_process_solutions[pcs_idx];
            _output->doOutputLastTimestep(pcs, spd->process_output, timestep, t,
                                          x);

            ++pcs_idx;
        }
    }

    return nonlinear_solver_succeeded;
}

bool UncoupledProcessesTimeLoop::solveUncoupledEquationSystems(
    const double t, const double dt, const std::size_t timestep_id)
{
    // TODO use process name
    unsigned pcs_idx = 0;
    for (auto& spd : _per_process_data)
    {
        auto& pcs = spd->process;
        BaseLib::RunTime time_timestep_process;
        time_timestep_process.start();

        auto& x = *_process_solutions[pcs_idx];
        pcs.preTimestep(x, t, dt);

        const auto nonlinear_solver_succeeded = solveOneTimeStepOneProcess(
            x, timestep_id, t, dt, *spd,
            ProcessLib::createVoidStaggeredCouplingTerm(), *_output);
        pcs.postTimestep(x);
        pcs.computeSecondaryVariable(t, x);

        INFO("[time] Solving process #%u took %g s in time step #%u ", pcs_idx,
             time_timestep_process.elapsed(), timestep_id);

        if (!nonlinear_solver_succeeded)
        {
            ERR("The nonlinear solver failed in time step #%u at t = %g "
                "s"
                " for process #%u.",
                timestep_id, t, pcs_idx);

            // save unsuccessful solution
            _output->doOutputAlways(pcs, spd->process_output, timestep_id, t,
                                    x);

            return false;
        }
        else
        {
            _output->doOutput(pcs, spd->process_output, timestep_id, t, x);
        }

        ++pcs_idx;
    }

    return true;
}

bool UncoupledProcessesTimeLoop::solveCoupledEquationSystemsByStaggeredScheme(
    const double t, const double dt, const std::size_t timestep_id)
{
    // Coupling iteration
    bool coupling_iteration_converged = true;
    for (unsigned i = 0; i < _global_coupling_max_iterations; i++)
    {
        // TODO use process name
        bool nonlinear_solver_succeeded = true;
        unsigned pcs_idx = 0;
        for (auto& spd : _per_process_data)
        {
            auto& pcs = spd->process;
            BaseLib::RunTime time_timestep_process;
            time_timestep_process.start();

            auto& x = *_process_solutions[pcs_idx];
            if (i == 0)
            {
                pcs.preTimestep(x, t, dt);

                // Create a coupling term
                _global_coupling_conv_crit->preFirstIteration();
            }
            else
            {
                _global_coupling_conv_crit->setNoFirstIteration();
            }
            StaggeredCouplingTerm coupled_term(
                spd->coupled_processes,
                _solutions_of_coupled_processes[pcs_idx], dt);

            const auto nonlinear_solver_succeeded = solveOneTimeStepOneProcess(
                x, timestep_id, t, dt, *spd, coupled_term, *_output);

            INFO(
                "[time] Solving process #%u took %g s in time step #%u "
                " coupling iteration #%u",
                pcs_idx, time_timestep_process.elapsed(), timestep_id, i);

            if (!nonlinear_solver_succeeded)
            {
                ERR("The nonlinear solver failed in time step #%u at t = %g "
                    "s"
                    " for process #%u.",
                    timestep_id, t, pcs_idx);

                // save unsuccessful solution
                _output->doOutputAlways(pcs, spd->process_output, timestep_id,
                                        t, x);

                break;
            }

            // Check the convergence of the coupling iteration
            auto& x_old = *_solutions_of_last_cpl_iteration[pcs_idx];
            MathLib::LinAlg::axpy(x_old, -1.0, x);
            _global_coupling_conv_crit->checkResidual(x_old);
            coupling_iteration_converged =
                coupling_iteration_converged &&
                _global_coupling_conv_crit->isSatisfied();

            if (coupling_iteration_converged)
                break;
            MathLib::LinAlg::copy(x, x_old);

            ++pcs_idx;
        }

        if (coupling_iteration_converged)
            break;

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
        auto& pcs = spd->process;
        auto& x = *_process_solutions[pcs_idx];
        pcs.postTimestep(x);
        pcs.computeSecondaryVariable(t, x);

        _output->doOutput(pcs, spd->process_output, timestep_id, t, x);
        ++pcs_idx;
    }

    return true;
}

UncoupledProcessesTimeLoop::~UncoupledProcessesTimeLoop()
{
    for (auto* x : _process_solutions)
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);

    for (auto* x : _solutions_of_last_cpl_iteration)
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
}

}  // namespace ProcessLib
