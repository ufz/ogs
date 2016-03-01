/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef APPLICATIONSLIB_UNCOUPLED_PROCESSES_TIMELOOP
#define APPLICATIONSLIB_UNCOUPLED_PROCESSES_TIMELOOP

#include <memory>

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "NumLib/ODESolver/NonlinearSolver.h"

#include "ProjectData.h"


namespace ApplicationsLib
{

//! Time loop capable of time-integrating several uncoupled processes at once.
template<typename Matrix, typename Vector>
class UncoupledProcessesTimeLoop
{
public:
    bool loop(ProjectData& project, std::string const& outdir);

private:
    //! An abstract nonlinear solver
    using AbstractNLSolver = NumLib::NonlinearSolverBase<Matrix, Vector>;
    //! An abstract equations system
    using EquationSystem   = NumLib::EquationSystem;
    //! An abstract process
    using Process          = ProcessLib::Process<GlobalSetupType>;
    //! An abstract time discretization
    using TimeDisc         = NumLib::TimeDiscretization<Vector>;

    struct SingleProcessData
    {
        template<NumLib::ODESystemTag ODETag, NumLib::NonlinearSolverTag NLTag>
        SingleProcessData(
                NumLib::NonlinearSolver<Matrix, Vector, NLTag>& nonlinear_solver,
                TimeDisc& time_disc,
                NumLib::ODESystem<Matrix, Vector, ODETag, NLTag>& ode_sys)
            : nonlinear_solver_tag(NLTag)
            , nonlinear_solver(nonlinear_solver)
            , tdisc_ode_sys(
                  new NumLib::TimeDiscretizedODESystem<Matrix, Vector, ODETag, NLTag>(
                                ode_sys, time_disc))
            , mat_strg(dynamic_cast<NumLib::InternalMatrixStorage&>(*tdisc_ode_sys))
        {}

        //! Tag containing the missing type information necessary to cast the other
        //! members of this struct to their concrety types.
        NumLib::NonlinearSolverTag const nonlinear_solver_tag;
        AbstractNLSolver& nonlinear_solver;
        //! type-erased time-discretized ODE system
        std::unique_ptr<EquationSystem> tdisc_ode_sys;
        //! cast of \c tdisc_ode_sys to NumLib::InternalMatrixStorage
        NumLib::InternalMatrixStorage& mat_strg;
    };

    /*! Creates a new instance of PerProcessData from the given data.
     *
     * \tparam ODETag tag indicating the type of ODE to be solved.
     *
     * \param nonlinear_solver the nonlinear solver to be used
     * \param ode_sys   the ODE (i.e. the process) to be integrated in time
     * \param time_disc the time discretization used for integrating \c ode_sys
     *
     * \note Currently the \c ODETag is automatically inferred from the given
     *       \c ody_sys. This works as long as \c Process derives from
     *       \c ODESystem<Matrix, Vector, ODETag>, i.e. as long we only deal with
     *       one type of ODE. When we introduce more types, this method will have
     *       to be extended slightly.
     */
    template<NumLib::ODESystemTag ODETag>
    static
    SingleProcessData
    makeSingleProcessData(
            AbstractNLSolver& nonlinear_solver,
            NumLib::ODESystem<Matrix, Vector, ODETag,
                NumLib::NonlinearSolverTag::Picard>& ode_sys,
            TimeDisc& time_disc)
    {
        using Tag = NumLib::NonlinearSolverTag;
        // A concrete Picard solver
        using NonlinearSolverPicard =
            NumLib::NonlinearSolver<Matrix, Vector, Tag::Picard>;
        // A concrete Newton solver
        using NonlinearSolverNewton =
            NumLib::NonlinearSolver<Matrix, Vector, Tag::Newton>;

        if (auto* nonlinear_solver_picard =
            dynamic_cast<NonlinearSolverPicard*>(&nonlinear_solver))
        {
            // The Picard solver can also work with a Newton-ready ODE,
            // because the Newton ODESystem derives from the Picard ODESystem.
            // So no further checks are needed here.

            return SingleProcessData{ *nonlinear_solver_picard, time_disc, ode_sys };
        }
        else if (auto* nonlinear_solver_newton =
                 dynamic_cast<NonlinearSolverNewton*>(&nonlinear_solver))
        {
            // The Newton-Raphson method needs a Newton-ready ODE.

            using ODENewton = NumLib::ODESystem<Matrix, Vector, ODETag, Tag::Newton>;
            if (auto* ode_newton = dynamic_cast<ODENewton*>(&ode_sys))
            {
                return SingleProcessData{
                    *nonlinear_solver_newton, time_disc, *ode_newton };
            }
            else {
                ERR("You are trying to solve a non-Newton-ready ODE with the"
                    " Newton-Raphson method. Aborting");
                std::abort();
            }
        }
        else {
            ERR("Encountered unknown nonlinear solver type. Aborting");
            std::abort();
        }
    }

    //! Constructs and returns a \c vector of PerProcessData from the given
    //! ProjectData.
    std::vector<SingleProcessData> initInternalData(ProjectData& project);

    //! Sets and returns initial conditions for the given ProjectData and
    //! PerProcessData.
    std::vector<Vector> setInitialConditions(
            ProjectData& project, double const t0,
            std::vector<SingleProcessData>& per_process_data);

    //! Solves one timestep for the given \c process.
    bool solveOneTimeStepOneProcess(
            unsigned const timestep,
            Vector& x, double const t, double const delta_t,
            SingleProcessData& process_data,
            Process& process,
            std::string const& output_file_name);

    //! Sets the EquationSystem for the given nonlinear solver,
    //! which is Picard or Newton depending on the NLTag.
    template<NumLib::NonlinearSolverTag NLTag>
    static void setEquationSystem(AbstractNLSolver& nonlinear_solver,
                                  EquationSystem& eq_sys)
    {
        using Solver = NumLib::NonlinearSolver<Matrix, Vector, NLTag>;
        using EqSys  = NumLib::NonlinearSystem<Matrix, Vector, NLTag>;

        assert(dynamic_cast<Solver*>(&nonlinear_solver) != nullptr);
        assert(dynamic_cast<EqSys*> (&eq_sys) != nullptr);

        auto& nl_solver_ = static_cast<Solver&>(nonlinear_solver);
        auto& eq_sys_    = static_cast<EqSys&> (eq_sys);

        nl_solver_.setEquationSystem(eq_sys_);
    }

    //! Sets the EquationSystem for the given nonlinear solver,
    //! transparently both for Picard and Newton solvers.
    static void setEquationSystem(
            AbstractNLSolver& nonlinear_solver, EquationSystem& eq_sys,
            NumLib::NonlinearSolverTag nl_tag)
    {
        using Tag = NumLib::NonlinearSolverTag;
        switch (nl_tag)
        {
        case Tag::Picard:
            setEquationSystem<Tag::Picard>(nonlinear_solver, eq_sys);
            break;
        case Tag::Newton:
            setEquationSystem<Tag::Newton>(nonlinear_solver, eq_sys);
            break;
        }
    }
};

//! Builds an UncoupledProcessesTimeLoop from the given configuration.
template<typename Matrix, typename Vector>
std::unique_ptr<UncoupledProcessesTimeLoop<Matrix, Vector> >
createUncoupledProcessesTimeLoop(BaseLib::ConfigTree const& conf)
{
    auto const type = conf.getConfParam<std::string>("type");

    if (type == "SingleStep")
    {
        using TimeLoop = UncoupledProcessesTimeLoop<Matrix, Vector>;
        return std::unique_ptr<TimeLoop>(new TimeLoop);
    }
    else
    {
            ERR("Unknown timestepper type: `%s'.", type.c_str());
            std::abort();
    }
}


template<typename Matrix, typename Vector>
std::vector<typename UncoupledProcessesTimeLoop<Matrix, Vector>::SingleProcessData>
UncoupledProcessesTimeLoop<Matrix, Vector>::
initInternalData(ProjectData& project)
{
    auto const num_processes = std::distance(project.processesBegin(),
                                             project.processesEnd());

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


template<typename Matrix, typename Vector>
std::vector<Vector>
UncoupledProcessesTimeLoop<Matrix, Vector>::
setInitialConditions(ProjectData& project,
                     double const t0,
                     std::vector<SingleProcessData>& per_process_data)
{
    auto const num_processes = std::distance(project.processesBegin(),
                                             project.processesEnd());

    std::vector<Vector> process_solutions;
    process_solutions.reserve(num_processes);

    unsigned pcs_idx = 0;
    for (auto p = project.processesBegin(); p != project.processesEnd();
         ++p, ++pcs_idx)
    {
        auto& pcs         = **p;
        auto& time_disc   =   pcs.getTimeDiscretization();

        auto& ppd         =   per_process_data[pcs_idx];
        auto& ode_sys     =  *ppd.tdisc_ode_sys;
        auto const nl_tag =   ppd.nonlinear_solver_tag;

        auto const num_eqs = ode_sys.getNumEquations();

        // TODO maybe more is required for PETSc
        // append a solution vector of size num_eqs
        process_solutions.emplace_back(num_eqs);

        auto& x0 = process_solutions[pcs_idx];
        pcs.setInitialConditions(x0);

        time_disc.setInitialState(t0, x0); // push IC

        if (time_disc.needsPreload())
        {
            auto& nonlinear_solver = ppd.nonlinear_solver;
            auto& mat_strg         = ppd.mat_strg;

            setEquationSystem(nonlinear_solver, ode_sys, nl_tag);
            nonlinear_solver.assemble(x0);
            time_disc.pushState(t0, x0, mat_strg); // TODO: that might do duplicate work
        }
    }

    return process_solutions;
}


template<typename Matrix, typename Vector>
bool
UncoupledProcessesTimeLoop<Matrix, Vector>::
solveOneTimeStepOneProcess(
        unsigned const timestep,
        Vector& x, double const t, double const delta_t,
        SingleProcessData& process_data,
        typename UncoupledProcessesTimeLoop<Matrix, Vector>::Process& process,
        std::string const& output_file_name)
{
    auto& time_disc        =  process.getTimeDiscretization();
    auto& ode_sys          = *process_data.tdisc_ode_sys;
    auto& nonlinear_solver =  process_data.nonlinear_solver;
    auto const nl_tag      =  process_data.nonlinear_solver_tag;

    setEquationSystem(nonlinear_solver, ode_sys, nl_tag);

    // INFO("time: %e, delta_t: %e", t, delta_t);
    time_disc.nextTimestep(t, delta_t);

    bool nonlinear_solver_succeeded = nonlinear_solver.solve(x);

    auto& mat_strg = process_data.mat_strg;
    time_disc.pushState(t, x, mat_strg);

    process.postTimestep(output_file_name, timestep, x);

    return nonlinear_solver_succeeded;
}


template<typename Matrix, typename Vector>
bool
UncoupledProcessesTimeLoop<Matrix, Vector>::
loop(ProjectData& project, std::string const& outdir)
{
    auto per_process_data = initInternalData(project);

    auto const t0 = 0.0; // time of the IC

    // init solution storage
    auto process_solutions = setInitialConditions(project, t0, per_process_data);

    auto const delta_t = 1.0;

    // TODO only for now
    // Make sure there will be exactly one iteration of the loop below.
    auto const t_end   = 1.5;

    double t;
    unsigned timestep = 1; // the first timestep really is number one
    bool nonlinear_solver_succeeded = true;
    for (t=t0+delta_t; t<t_end+delta_t; t+=delta_t, ++timestep)
    {
        unsigned pcs_idx = 0;
        for (auto p = project.processesBegin(); p != project.processesEnd();
             ++p, ++pcs_idx)
        {
            std::string const& outpref = project.getOutputFilePrefix();
            std::string const  output_file_name =
                    BaseLib::joinPaths(outdir, outpref)
                    + "_pcs_" + std::to_string(pcs_idx)
                    + "_ts_"  + std::to_string(timestep)
                    // + "_t_"   + std::to_string(t) // TODO: add that later
                    + ".vtu";

            auto& x = process_solutions[pcs_idx];

            nonlinear_solver_succeeded = solveOneTimeStepOneProcess(
                        timestep, x, t, delta_t,
                        per_process_data[pcs_idx], **p, output_file_name);

            if (!nonlinear_solver_succeeded) {
                ERR("The nonlinear solver failed in timestep #%u at t = %g s"
                    " for process #%u.", timestep, t, pcs_idx);
                break;
            }
        }

        if (!nonlinear_solver_succeeded) break;

        break; // TODO only do a single timestep for now
    }

    return nonlinear_solver_succeeded;
}

}

#endif // APPLICATIONSLIB_UNCOUPLED_PROCESSES_TIMELOOP
