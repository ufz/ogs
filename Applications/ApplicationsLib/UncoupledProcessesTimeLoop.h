/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// TODO
#pragma once

#include <memory>

#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "NumLib/ODESolver/NonlinearSolver.h"

#include "ProjectData.h"

#include "BaseLib/ConfigTree.h"


namespace ApplicationsLib
{

template<typename Matrix, typename Vector, NumLib::NonlinearSolverTag NLTag>
class UncoupledProcessesTimeLoop
{
public:
    using TDiscODESys = NumLib::TimeDiscretizedODESystemBase<Matrix, Vector, NLTag>;
    using NLSolver = NumLib::NonlinearSolver<Matrix, Vector, NLTag>;

    /*! Constructs an new instance.
     * \param ode_sys The ODE system to be integrated
     * \param nonlinear_solver The solver to be used to resolve nonlinearities.
     */
    explicit
    UncoupledProcessesTimeLoop(
            std::unique_ptr<NLSolver>&& nonlinear_solver
            )
        : _nonlinear_solver(std::move(nonlinear_solver))
    {}

    bool loop(ProjectData& project);

private:
    using TimeDisc = NumLib::TimeDiscretization<Vector>;

    //! \note both members of this struct are pointers to abstract types.
    struct PerProcessData
    {
        std::unique_ptr<TimeDisc>    time_disc;
        std::unique_ptr<TDiscODESys> tdisc_ode_sys;
    };


    template<NumLib::ODESystemTag ODETag>
    static
    PerProcessData
    makePerProcessData(
            NumLib::ODESystem<Matrix, Vector, ODETag, NLTag>& ode_sys,
            std::unique_ptr<TimeDisc>&& time_disc)
    {
        PerProcessData ppd{ std::move(time_disc) };
        ppd.tdisc_ode_sys.reset(
            new NumLib::TimeDiscretizedODESystem<Matrix, Vector, ODETag, NLTag>(
                ode_sys, *ppd.time_disc));

        return ppd;
    }


    std::unique_ptr<NLSolver> _nonlinear_solver;
};

template<typename Matrix, typename Vector, NumLib::NonlinearSolverTag NLTag>
std::unique_ptr<UncoupledProcessesTimeLoop<Matrix, Vector, NLTag> >
createUncoupledProcessesTimeLoop(BaseLib::ConfigTree const& conf)
{
    auto const type = conf.getConfParam<std::string>("type");

    if (type == "SingleStep")
    {
        using TimeLoop = UncoupledProcessesTimeLoop<Matrix, Vector, NLTag>;
        using NLSolver = typename TimeLoop::NLSolver;

        auto const tol = 1e-6;
        auto const max_iter = 10u;

        std::unique_ptr<NLSolver> nl_solver(new NLSolver(tol, max_iter));

        return std::unique_ptr<TimeLoop>(new TimeLoop(std::move(nl_solver)));
    }
    else
    {
            ERR("Unknown timestepper type: `%s'.", type.c_str());
            std::abort();
    }
}




template<typename Matrix, typename Vector, NumLib::NonlinearSolverTag NLTag>
bool
UncoupledProcessesTimeLoop<Matrix, Vector, NLTag>::
loop(ProjectData& project)
{

#if 0
    unsigned i = 0;  // process counter, used to distinguish output files
    for (auto p = project.processesBegin(); p != project.processesEnd();
         ++p)
    {
        accepted = accepted && (*p)->solve(dt);

        if (!accepted)
        {
            ERR("Timestep has not been accepted. Aborting.");
            break;
        }

        std::string const output_file_name =
                BaseLib::joinPaths(outdir, out_pref) +
                "_pcs_" + std::to_string(i) + "_ts_" +
                std::to_string(timestep) + ".vtu";
        (*p)->postTimestep(output_file_name, timestep);

        ++i;
    }

#endif

    auto const num_processes = std::distance(project.processesBegin(),
                                             project.processesEnd());

    std::vector<PerProcessData> time_disc_ode_syss;
    time_disc_ode_syss.reserve(num_processes);

    // create a time discretized ODE system for each process
    for (auto p = project.processesBegin(); p != project.processesEnd();
         ++p)
    {
        auto time_disc = std::unique_ptr<TimeDisc>(
                            new NumLib::BackwardEuler<Vector>);
        time_disc_ode_syss.emplace_back(
                    makePerProcessData(**p, std::move(time_disc)));
    }

    // init solution storage
    std::vector<Vector> process_solutions(num_processes); // TODO: waste of memory

    // set ICs
#if 0
    auto& time_disc = _ode_sys.getTimeDiscretization();

    time_disc.setInitialState(t0, x0); // push IC

    if (time_disc.needsPreload()) {
        _nonlinear_solver.assemble(_ode_sys, x);
        time_disc.pushState(t0, x0, _ode_sys); // TODO: that might do duplicate work
    }
#endif
    // Vector x(x0); // solution vector

    auto const t0      = 0.0;
    auto const delta_t = 1.0;

    // Make sure there will be exactly one iteration of the loop below.
    auto const t_end   = 1.5;

    double t;
    unsigned timestep = 0;
    bool nl_slv_succeeded = true;
    for (t=t0+delta_t; t<t_end+delta_t; t+=delta_t, ++timestep)
    {

        unsigned pcs_idx = 0;
        for (auto p = project.processesBegin(); p != project.processesEnd();
             ++p, ++pcs_idx)
        {
            auto& time_disc = *time_disc_ode_syss[pcs_idx].time_disc;
            auto& ode_sys   = *time_disc_ode_syss[pcs_idx].tdisc_ode_sys;
            auto& x = process_solutions[pcs_idx];

            // INFO("time: %e, delta_t: %e", t, delta_t);
            time_disc.nextTimestep(t, delta_t);

            nl_slv_succeeded = _nonlinear_solver->solve(ode_sys, x);

            // TODO error message
            if (!nl_slv_succeeded) break;

            time_disc.pushState(t, x, ode_sys);

            auto const  t_cb = t; // make sure the callback cannot overwrite anything.
            auto const& x_cb = x; // ditto.

            // TODO
            // post_timestep(t_cb, x_cb);
        }


        break; // TODO only do a single timestep for now
    }

    if (!nl_slv_succeeded) {
        ERR("Nonlinear solver failed in timestep #%u at t = %g s", timestep, t);
    }
    return nl_slv_succeeded;
}

}
