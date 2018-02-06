/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/ODESolver/Types.h"
#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"

#include "ProcessLib/Output/ProcessOutput.h"

#include "CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
class Process;
}

namespace ProcessLib
{
struct ProcessData
{
    template <NumLib::NonlinearSolverTag NLTag>
    ProcessData(std::unique_ptr<NumLib::TimeStepAlgorithm>&& timestepper_,
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

    ProcessData(ProcessData&& pd)
        : timestepper(std::move(pd.timestepper)),
          nonlinear_solver_tag(pd.nonlinear_solver_tag),
          nonlinear_solver(pd.nonlinear_solver),
          nonlinear_solver_converged(pd.nonlinear_solver_converged),
          conv_crit(std::move(pd.conv_crit)),
          time_disc(std::move(pd.time_disc)),
          tdisc_ode_sys(std::move(pd.tdisc_ode_sys)),
          mat_strg(pd.mat_strg),
          process(pd.process),
          process_output(std::move(pd.process_output))
    {
        pd.mat_strg = nullptr;
    }

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

    /// Process ID. It is alway 0 when the monolithic scheme is used or
    /// a single process is modelled.
    int process_id = 0;

    Process& process;
    ProcessOutput process_output;
};
}  // namespace ProcessLib
