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
#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

#include "ProjectData.h"


namespace ApplicationsLib
{

//! Time loop capable of time-integrating several uncoupled processes at once.
class UncoupledProcessesTimeLoop
{
public:
    explicit UncoupledProcessesTimeLoop(
            std::unique_ptr<NumLib::ITimeStepAlgorithm>&& timestepper)
        : _timestepper{std::move(timestepper)}
    {}

    bool loop(ProjectData& project);

    ~UncoupledProcessesTimeLoop();

private:
    //! An abstract nonlinear solver
    using AbstractNLSolver = NumLib::NonlinearSolverBase<GlobalMatrix, GlobalVector>;
    //! An abstract equations system
    using EquationSystem   = NumLib::EquationSystem<GlobalVector>;
    //! An abstract process
    using Process          = ProcessLib::Process;
    //! An abstract time discretization
    using TimeDisc         = NumLib::TimeDiscretization<GlobalVector>;

    std::vector<GlobalVector*> _process_solutions;
    std::unique_ptr<NumLib::ITimeStepAlgorithm> _timestepper;

    struct SingleProcessData
    {
        template<NumLib::ODESystemTag ODETag, NumLib::NonlinearSolverTag NLTag>
        SingleProcessData(
                NumLib::NonlinearSolver<GlobalMatrix, GlobalVector, NLTag>& nonlinear_solver,
                TimeDisc& time_disc,
                NumLib::ODESystem<GlobalMatrix, GlobalVector, ODETag, NLTag>& ode_sys)
            : nonlinear_solver_tag(NLTag)
            , nonlinear_solver(nonlinear_solver)
            , tdisc_ode_sys(
                  new NumLib::TimeDiscretizedODESystem<GlobalMatrix, GlobalVector, ODETag, NLTag>(
                                ode_sys, time_disc))
            , mat_strg(dynamic_cast<NumLib::InternalMatrixStorage&>(*tdisc_ode_sys))
        {}

        SingleProcessData(SingleProcessData&& spd)
            : nonlinear_solver_tag(spd.nonlinear_solver_tag)
            , nonlinear_solver    (spd.nonlinear_solver)
            , tdisc_ode_sys       (std::move(spd.tdisc_ode_sys))
            , mat_strg            (spd.mat_strg)
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
     *       \c ODESystem<Matrix, GlobalVector, ODETag>, i.e. as long we only deal with
     *       one type of ODE. When we introduce more types, this method will have
     *       to be extended slightly.
     */
    template<NumLib::ODESystemTag ODETag>
    SingleProcessData
    makeSingleProcessData(
            AbstractNLSolver& nonlinear_solver,
            NumLib::ODESystem<GlobalMatrix, GlobalVector, ODETag,
                NumLib::NonlinearSolverTag::Picard>& ode_sys,
            TimeDisc& time_disc)
    {
        using Tag = NumLib::NonlinearSolverTag;
        // A concrete Picard solver
        using NonlinearSolverPicard =
            NumLib::NonlinearSolver<GlobalMatrix, GlobalVector, Tag::Picard>;
        // A concrete Newton solver
        using NonlinearSolverNewton =
            NumLib::NonlinearSolver<GlobalMatrix, GlobalVector, Tag::Newton>;

        if (auto* nonlinear_solver_picard =
            dynamic_cast<NonlinearSolverPicard*>(&nonlinear_solver))
        {
            // The Picard solver can also work with a Newton-ready ODE,
            // because the Newton ODESystem derives from the Picard ODESystem.
            // So no further checks are needed here.

            return SingleProcessData{
                *nonlinear_solver_picard, time_disc, ode_sys };
        }
        else if (auto* nonlinear_solver_newton =
                 dynamic_cast<NonlinearSolverNewton*>(&nonlinear_solver))
        {
            // The Newton-Raphson method needs a Newton-ready ODE.

            using ODENewton = NumLib::ODESystem<GlobalMatrix, GlobalVector, ODETag, Tag::Newton>;
            if (auto* ode_newton = dynamic_cast<ODENewton*>(&ode_sys))
            {
                return SingleProcessData{
                    *nonlinear_solver_newton, time_disc, *ode_newton };
            }
            else {
                OGS_FATAL("You are trying to solve a non-Newton-ready ODE with the"
                    " Newton-Raphson method. Aborting");
            }
        }
        else {
            OGS_FATAL("Encountered unknown nonlinear solver type. Aborting");
        }
    }

    //! Constructs and returns a \c vector of SingleProcessData from the given
    //! \c project.
    std::vector<SingleProcessData> initInternalData(ProjectData& project);

    //! Sets initial conditions for the given \c project and \c per_process_data.
    void setInitialConditions(
            ProjectData& project, double const t0,
            std::vector<SingleProcessData>& per_process_data);

    //! Solves one timestep for the given \c process.
    bool solveOneTimeStepOneProcess(
            GlobalVector& x, double const t, double const delta_t,
            SingleProcessData& process_data,
            Process& process);

    //! Sets the EquationSystem for the given nonlinear solver,
    //! which is Picard or Newton depending on the NLTag.
    template<NumLib::NonlinearSolverTag NLTag>
    static void setEquationSystem(AbstractNLSolver& nonlinear_solver,
                                  EquationSystem& eq_sys)
    {
        using Solver = NumLib::NonlinearSolver<GlobalMatrix, GlobalVector, NLTag>;
        using EqSys  = NumLib::NonlinearSystem<GlobalMatrix, GlobalVector, NLTag>;

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

    //! Applies known solutions to the solution vector \c x, transparently
    //! for equation systems linearized with either the Picard or Newton method.
    template<NumLib::NonlinearSolverTag NLTag>
    static void applyKnownSolutions(
            EquationSystem const& eq_sys, GlobalVector& x)
    {
        using EqSys = NumLib::NonlinearSystem<GlobalMatrix, GlobalVector, NLTag>;
        assert(dynamic_cast<EqSys const*> (&eq_sys) != nullptr);
        auto& eq_sys_ = static_cast<EqSys const&> (eq_sys);

        eq_sys_.applyKnownSolutions(x);
    }

    //! Applies known solutions to the solution vector \c x, transparently
    //! for equation systems linearized with either the Picard or Newton method.
    static void applyKnownSolutions(
            EquationSystem const& eq_sys,
            NumLib::NonlinearSolverTag const nl_tag, GlobalVector& x)
    {
        using Tag = NumLib::NonlinearSolverTag;
        switch (nl_tag)
        {
        case Tag::Picard: applyKnownSolutions<Tag::Picard>(eq_sys, x); break;
        case Tag::Newton: applyKnownSolutions<Tag::Newton>(eq_sys, x); break;
        }
    }
};

//! Builds an UncoupledProcessesTimeLoop from the given configuration.
inline std::unique_ptr<UncoupledProcessesTimeLoop>
createUncoupledProcessesTimeLoop(BaseLib::ConfigTree const& conf)
{
    //! \ogs_file_param{prj__time_stepping__type}
    auto const type = conf.peekConfigParameter<std::string>("type");

    std::unique_ptr<NumLib::ITimeStepAlgorithm> timestepper;

    if (type == "SingleStep") {
        conf.ignoreConfigParameter("type");
        timestepper.reset(new NumLib::FixedTimeStepping(0.0, 1.0, 1.0));
    } else if (type == "FixedTimeStepping") {
        timestepper = NumLib::FixedTimeStepping::newInstance(conf);
    } else {
            OGS_FATAL("Unknown timestepper type: `%s'.", type.c_str());
    }

    using TimeLoop = UncoupledProcessesTimeLoop;
    return std::unique_ptr<TimeLoop>{new TimeLoop{std::move(timestepper)}};
}


inline std::vector<typename UncoupledProcessesTimeLoop::SingleProcessData>
UncoupledProcessesTimeLoop::
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


inline void
UncoupledProcessesTimeLoop::
setInitialConditions(ProjectData& project,
                     double const t0,
                     std::vector<SingleProcessData>& per_process_data)
{
    auto const num_processes = std::distance(project.processesBegin(),
                                             project.processesEnd());

    _process_solutions.reserve(num_processes);

    unsigned pcs_idx = 0;
    for (auto p = project.processesBegin(); p != project.processesEnd();
         ++p, ++pcs_idx)
    {
        auto& pcs         = **p;
        auto& time_disc   =   pcs.getTimeDiscretization();

        auto& ppd         =   per_process_data[pcs_idx];
        auto& ode_sys     =  *ppd.tdisc_ode_sys;
        auto const nl_tag =   ppd.nonlinear_solver_tag;

        // append a solution vector of suitable size
        _process_solutions.emplace_back(
                    &MathLib::GlobalVectorProvider<GlobalVector>::provider.getVector(
                        ode_sys.getMatrixSpecifications()));

        auto& x0 = *_process_solutions[pcs_idx];
        pcs.setInitialConditions(x0);
        MathLib::BLAS::finalizeAssembly(x0);

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
}


inline bool
UncoupledProcessesTimeLoop::
solveOneTimeStepOneProcess(
        GlobalVector& x, double const t, double const delta_t,
        SingleProcessData& process_data,
        typename UncoupledProcessesTimeLoop::Process& process)
{
    auto& time_disc        =  process.getTimeDiscretization();
    auto& ode_sys          = *process_data.tdisc_ode_sys;
    auto& nonlinear_solver =  process_data.nonlinear_solver;
    auto const nl_tag      =  process_data.nonlinear_solver_tag;

    setEquationSystem(nonlinear_solver, ode_sys, nl_tag);

    // Note: Order matters!
    // First advance to the next timestep, then set known solutions at that
    // time, afterwards pass the right solution vector and time to the
    // preTimestep() hook.

    time_disc.nextTimestep(t, delta_t);

    applyKnownSolutions(ode_sys, nl_tag, x);

    process.preTimestep(x, t, delta_t);

    bool nonlinear_solver_succeeded = nonlinear_solver.solve(x);

    auto& mat_strg = process_data.mat_strg;
    time_disc.pushState(t, x, mat_strg);

    process.postTimestep(x);

    return nonlinear_solver_succeeded;
}


inline bool
UncoupledProcessesTimeLoop::
loop(ProjectData& project)
{
    auto per_process_data = initInternalData(project);

    auto& out_ctrl = project.getOutputControl();
    out_ctrl.initialize(project.processesBegin(), project.processesEnd());

    auto const t0 = _timestepper->getTimeStep().current(); // time of the IC

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
    std::size_t timestep = 1; // the first timestep really is number one
    bool nonlinear_solver_succeeded = true;

    while (_timestepper->next())
    {
        auto const ts = _timestepper->getTimeStep();
        auto const delta_t = ts.dt();
        t        = ts.current();
        timestep = ts.steps();

        INFO("=== timestep #%u (t=%gs, dt=%gs) ==============================",
             timestep, t, delta_t);

        // TODO use process name
        unsigned pcs_idx = 0;
        for (auto p = project.processesBegin(); p != project.processesEnd();
             ++p, ++pcs_idx)
        {
            auto& x = *_process_solutions[pcs_idx];

            nonlinear_solver_succeeded = solveOneTimeStepOneProcess(
                        x, t, delta_t, per_process_data[pcs_idx], **p);

            if (!nonlinear_solver_succeeded) {
                ERR("The nonlinear solver failed in timestep #%u at t = %g s"
                    " for process #%u.", timestep, t, pcs_idx);

                // save unsuccessful solution
                out_ctrl.doOutputAlways(**p, timestep, t, x);

                break;
            } else {
                out_ctrl.doOutput(**p, timestep, t, x);
            }
        }

        if (!nonlinear_solver_succeeded) break;
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


inline UncoupledProcessesTimeLoop::
~UncoupledProcessesTimeLoop()
{
    for (auto * x : _process_solutions)
        MathLib::GlobalVectorProvider<GlobalVector>::provider.releaseVector(*x);
}

} // namespace ApplicationsLib

#endif // APPLICATIONSLIB_UNCOUPLED_PROCESSES_TIMELOOP
