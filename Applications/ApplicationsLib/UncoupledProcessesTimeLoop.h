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
#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
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
    {
    }

    bool loop(ProjectData& project);

    ~UncoupledProcessesTimeLoop();

private:
    //! An abstract nonlinear solver
    using AbstractNLSolver = NumLib::NonlinearSolverBase;
    //! An abstract equations system
    using EquationSystem = NumLib::EquationSystem;
    //! An abstract process
    using Process = ProcessLib::Process;
    //! An abstract time discretization
    using TimeDisc = NumLib::TimeDiscretization;

    std::vector<GlobalVector*> _process_solutions;
    std::unique_ptr<NumLib::ITimeStepAlgorithm> _timestepper;

    struct SingleProcessData
    {
        template <NumLib::ODESystemTag ODETag, NumLib::NonlinearSolverTag NLTag>
        SingleProcessData(NumLib::NonlinearSolver<NLTag>& nonlinear_solver,
                          TimeDisc& time_disc,
                          NumLib::ODESystem<ODETag, NLTag>& ode_sys)
            : nonlinear_solver_tag(NLTag),
              nonlinear_solver(nonlinear_solver),
              tdisc_ode_sys(new NumLib::TimeDiscretizedODESystem<ODETag, NLTag>(
                  ode_sys, time_disc)),
              mat_strg(
                  dynamic_cast<NumLib::InternalMatrixStorage&>(*tdisc_ode_sys))
        {
        }

        SingleProcessData(SingleProcessData&& spd)
            : nonlinear_solver_tag(spd.nonlinear_solver_tag),
              nonlinear_solver(spd.nonlinear_solver),
              tdisc_ode_sys(std::move(spd.tdisc_ode_sys)),
              mat_strg(spd.mat_strg)
        {
        }

        //! Tag containing the missing type information necessary to cast the
        //! other members of this struct to their concrety types.
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
     *       \c ODESystem<GlobalMatrix, GlobalVector, ODETag>, i.e. as long we
     *       only deal with one type of ODE. When we introduce more types, this
     *       method will have to be extended slightly.
     */
    template <NumLib::ODESystemTag ODETag>
    SingleProcessData makeSingleProcessData(
        AbstractNLSolver& nonlinear_solver,
        NumLib::ODESystem<ODETag, NumLib::NonlinearSolverTag::Picard>& ode_sys,
        TimeDisc& time_disc)
    {
        using Tag = NumLib::NonlinearSolverTag;
        // A concrete Picard solver
        using NonlinearSolverPicard = NumLib::NonlinearSolver<Tag::Picard>;
        // A concrete Newton solver
        using NonlinearSolverNewton = NumLib::NonlinearSolver<Tag::Newton>;

        if (auto* nonlinear_solver_picard =
                dynamic_cast<NonlinearSolverPicard*>(&nonlinear_solver))
        {
            // The Picard solver can also work with a Newton-ready ODE,
            // because the Newton ODESystem derives from the Picard ODESystem.
            // So no further checks are needed here.

            return SingleProcessData{*nonlinear_solver_picard, time_disc,
                                     ode_sys};
        }
        else if (auto* nonlinear_solver_newton =
                     dynamic_cast<NonlinearSolverNewton*>(&nonlinear_solver))
        {
            // The Newton-Raphson method needs a Newton-ready ODE.

            using ODENewton = NumLib::ODESystem<ODETag, Tag::Newton>;
            if (auto* ode_newton = dynamic_cast<ODENewton*>(&ode_sys))
            {
                return SingleProcessData{*nonlinear_solver_newton, time_disc,
                                         *ode_newton};
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

    //! Constructs and returns a \c vector of SingleProcessData from the given
    //! \c project.
    std::vector<SingleProcessData> initInternalData(ProjectData& project);

    //! Sets initial conditions for the given \c project and \c
    //! per_process_data.
    void setInitialConditions(ProjectData& project, double const t0,
                              std::vector<SingleProcessData>& per_process_data);

    //! Solves one timestep for the given \c process.
    bool solveOneTimeStepOneProcess(
            GlobalVector& x, std::size_t const timestep, double const t, double const delta_t,
            SingleProcessData& process_data,
            Process& process, ProcessLib::Output const& output_control);

    //! Sets the EquationSystem for the given nonlinear solver,
    //! which is Picard or Newton depending on the NLTag.
    template <NumLib::NonlinearSolverTag NLTag>
    static void setEquationSystem(AbstractNLSolver& nonlinear_solver,
                                  EquationSystem& eq_sys,
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
    static void setEquationSystem(AbstractNLSolver& nonlinear_solver,
                                  EquationSystem& eq_sys,
                                  NumLib::ConvergenceCriterion& conv_crit,
                                  NumLib::NonlinearSolverTag nl_tag)
    {
        using Tag = NumLib::NonlinearSolverTag;
        switch (nl_tag)
        {
            case Tag::Picard:
                setEquationSystem<Tag::Picard>(nonlinear_solver, eq_sys,
                                               conv_crit);
                break;
            case Tag::Newton:
                setEquationSystem<Tag::Newton>(nonlinear_solver, eq_sys,
                                               conv_crit);
                break;
        }
    }

    //! Applies known solutions to the solution vector \c x, transparently
    //! for equation systems linearized with either the Picard or Newton method.
    template <NumLib::NonlinearSolverTag NLTag>
    static void applyKnownSolutions(EquationSystem const& eq_sys,
                                    GlobalVector& x)
    {
        using EqSys = NumLib::NonlinearSystem<NLTag>;
        assert(dynamic_cast<EqSys const*>(&eq_sys) != nullptr);
        auto& eq_sys_ = static_cast<EqSys const&>(eq_sys);

        eq_sys_.applyKnownSolutions(x);
    }

    //! Applies known solutions to the solution vector \c x, transparently
    //! for equation systems linearized with either the Picard or Newton method.
    static void applyKnownSolutions(EquationSystem const& eq_sys,
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
};

//! Builds an UncoupledProcessesTimeLoop from the given configuration.
std::unique_ptr<UncoupledProcessesTimeLoop> createUncoupledProcessesTimeLoop(
    BaseLib::ConfigTree const& conf);

}  // namespace ApplicationsLib

#endif  // APPLICATIONSLIB_UNCOUPLED_PROCESSES_TIMELOOP
