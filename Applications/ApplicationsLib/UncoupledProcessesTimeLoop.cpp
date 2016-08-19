/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UncoupledProcessesTimeLoop.h"
#include "BaseLib/RunTime.h"

namespace ApplicationsLib
{
std::unique_ptr<UncoupledProcessesTimeLoop> createUncoupledProcessesTimeLoop(
    BaseLib::ConfigTree const& conf)
{
    //! \ogs_file_param{prj__time_stepping__type}
    auto const type = conf.peekConfigParameter<std::string>("type");

    std::unique_ptr<NumLib::ITimeStepAlgorithm> timestepper;

    if (type == "SingleStep")
    {
        //! \ogs_file_param_special{prj__time_stepping__SingleStep}
        conf.ignoreConfigParameter("type");
        timestepper.reset(new NumLib::FixedTimeStepping(0.0, 1.0, 1.0));
    }
    else if (type == "FixedTimeStepping")
    {
        timestepper = NumLib::FixedTimeStepping::newInstance(conf);
    }
    else
    {
        OGS_FATAL("Unknown timestepper type: `%s'.", type.c_str());
    }

    using TimeLoop = UncoupledProcessesTimeLoop;
    return std::unique_ptr<TimeLoop>{new TimeLoop{std::move(timestepper)}};
}

std::vector<typename UncoupledProcessesTimeLoop::SingleProcessData>
UncoupledProcessesTimeLoop::initInternalData(ProjectData& project)
{
    auto const num_processes =
        std::distance(project.processesBegin(), project.processesEnd());

    std::vector<SingleProcessData> per_process_data;
    per_process_data.reserve(num_processes);

    // create a time discretized ODE system for each process
    for (auto p = project.processesBegin(); p != project.processesEnd(); ++p)
    {
        auto& pcs = **p;
        auto& nonlinear_solver = pcs.getNonlinearSolver();
        auto& time_disc = pcs.getTimeDiscretization();

        per_process_data.emplace_back(
            makeSingleProcessData(nonlinear_solver, **p, time_disc));
    }

    return per_process_data;
}

void UncoupledProcessesTimeLoop::setInitialConditions(
    ProjectData& project,
    double const t0,
    std::vector<SingleProcessData>& per_process_data)
{
    auto const num_processes =
        std::distance(project.processesBegin(), project.processesEnd());

    _process_solutions.reserve(num_processes);

    unsigned pcs_idx = 0;
    for (auto p = project.processesBegin(); p != project.processesEnd();
         ++p, ++pcs_idx)
    {
        auto& pcs = **p;
        auto& time_disc = pcs.getTimeDiscretization();

        auto& ppd = per_process_data[pcs_idx];
        auto& ode_sys = *ppd.tdisc_ode_sys;
        auto const nl_tag = ppd.nonlinear_solver_tag;

        // append a solution vector of suitable size
        _process_solutions.emplace_back(
            &NumLib::GlobalVectorProvider::provider.getVector(
                ode_sys.getMatrixSpecifications()));

        auto& x0 = *_process_solutions[pcs_idx];
        pcs.setInitialConditions(x0);
        MathLib::LinAlg::finalizeAssembly(x0);

        time_disc.setInitialState(t0, x0);  // push IC

        if (time_disc.needsPreload())
        {
            auto& nonlinear_solver = ppd.nonlinear_solver;
            auto& mat_strg = ppd.mat_strg;
            auto& conv_crit = pcs.getConvergenceCriterion();

            setEquationSystem(nonlinear_solver, ode_sys, conv_crit, nl_tag);
            nonlinear_solver.assemble(x0);
            time_disc.pushState(
                t0, x0, mat_strg);  // TODO: that might do duplicate work
        }
    }
}


bool
UncoupledProcessesTimeLoop::
solveOneTimeStepOneProcess(
        GlobalVector& x, std::size_t const timestep, double const t, double const delta_t,
        SingleProcessData& process_data,
        UncoupledProcessesTimeLoop::Process& process,
        ProcessLib::Output const& output_control)
{
    auto& time_disc = process.getTimeDiscretization();
    auto& conv_crit = process.getConvergenceCriterion();
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

    process.preTimestep(x, t, delta_t);

    auto const post_iteration_callback = [&](unsigned iteration,
                                             GlobalVector const& x) {
        output_control.doOutputNonlinearIteration(process, timestep, t, x, iteration);
    };

    bool nonlinear_solver_succeeded =
        nonlinear_solver.solve(x, post_iteration_callback);

    auto& mat_strg = process_data.mat_strg;
    time_disc.pushState(t, x, mat_strg);

    process.postTimestep(x);

    return nonlinear_solver_succeeded;
}

bool UncoupledProcessesTimeLoop::loop(ProjectData& project)
{
    auto per_process_data = initInternalData(project);

    auto& out_ctrl = project.getOutputControl();
    out_ctrl.initialize(project.processesBegin(), project.processesEnd());

    auto const t0 = _timestepper->getTimeStep().current();  // time of the IC

    // init solution storage
    setInitialConditions(project, t0, per_process_data);

    // output initial conditions
    {
        unsigned pcs_idx = 0;
        for (auto p = project.processesBegin(); p != project.processesEnd();
             ++p, ++pcs_idx)
        {
            auto const& x0 = *_process_solutions[pcs_idx];
            out_ctrl.doOutput(**p, 0, t0, x0);
        }
    }

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

        // TODO use process name
        unsigned pcs_idx = 0;
        for (auto p = project.processesBegin(); p != project.processesEnd();
             ++p, ++pcs_idx)
        {
            BaseLib::RunTime time_timestep_process;
            time_timestep_process.start();

            auto& x = *_process_solutions[pcs_idx];

            nonlinear_solver_succeeded = solveOneTimeStepOneProcess(
                x, timestep, t, delta_t, per_process_data[pcs_idx], **p,
                out_ctrl);

            INFO("[time] Solving process #%u took %g s in timestep #%u.",
                 timestep, time_timestep.elapsed(), pcs_idx);

            if (!nonlinear_solver_succeeded)
            {
                ERR("The nonlinear solver failed in timestep #%u at t = %g s"
                    " for process #%u.",
                    timestep, t, pcs_idx);

                // save unsuccessful solution
                out_ctrl.doOutputAlways(**p, timestep, t, x);

                break;
            }
            else
            {
                out_ctrl.doOutput(**p, timestep, t, x);
            }
        }

        INFO("[time] Timestep #%u took %g s.", timestep,
             time_timestep.elapsed());

        if (!nonlinear_solver_succeeded)
            break;
    }

    // output last timestep
    if (nonlinear_solver_succeeded)
    {
        unsigned pcs_idx = 0;
        for (auto p = project.processesBegin(); p != project.processesEnd();
             ++p, ++pcs_idx)
        {
            auto const& x = *_process_solutions[pcs_idx];
            out_ctrl.doOutputLastTimestep(**p, timestep, t, x);
        }
    }

    return nonlinear_solver_succeeded;
}

UncoupledProcessesTimeLoop::~UncoupledProcessesTimeLoop()
{
    for (auto* x : _process_solutions)
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
}

}  // namespace ApplicationsLib
