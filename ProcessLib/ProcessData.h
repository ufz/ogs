/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "CoupledSolutionsForStaggeredScheme.h"
#include "Process.h"

namespace ProcessLib
{
struct ProcessData
{
    ProcessData(std::unique_ptr<NumLib::TimeStepAlgorithm>&& timestepper_,
                NumLib::NonlinearSolverTag const nonlinear_solver_tag_,
                NumLib::NonlinearSolverBase& nonlinear_solver_,
                std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit_,
                std::unique_ptr<NumLib::TimeDiscretization>&& time_disc_,
                int const process_id_,
                Process& process_)
        : timestepper(std::move(timestepper_)),
          nonlinear_solver_tag(nonlinear_solver_tag_),
          nonlinear_solver(nonlinear_solver_),
          nonlinear_solver_status{true, 0},
          conv_crit(std::move(conv_crit_)),
          time_disc(std::move(time_disc_)),
          process_id(process_id_),
          process(process_)
    {
    }

    ProcessData(ProcessData&& pd)
        : timestepper(std::move(pd.timestepper)),
          nonlinear_solver_tag(pd.nonlinear_solver_tag),
          nonlinear_solver(pd.nonlinear_solver),
          nonlinear_solver_status(pd.nonlinear_solver_status),
          conv_crit(std::move(pd.conv_crit)),
          time_disc(std::move(pd.time_disc)),
          tdisc_ode_sys(std::move(pd.tdisc_ode_sys)),
          process_id(pd.process_id),
          process(pd.process)
    {
    }

    std::unique_ptr<NumLib::TimeStepAlgorithm> timestepper;

    //! Tag containing the missing type information necessary to cast the
    //! other members of this struct to their concrety types.
    NumLib::NonlinearSolverTag const nonlinear_solver_tag;
    NumLib::NonlinearSolverBase& nonlinear_solver;
    NumLib::NonlinearSolverStatus nonlinear_solver_status;
    std::unique_ptr<NumLib::ConvergenceCriterion> conv_crit;

    std::unique_ptr<NumLib::TimeDiscretization> time_disc;
    //! type-erased time-discretized ODE system
    std::unique_ptr<NumLib::EquationSystem> tdisc_ode_sys;

    int const process_id;

    Process& process;
};
}  // namespace ProcessLib
