/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/ODESolver/Types.h"
#include "NumLib/TimeStepping/Algorithms/TimeStepAlgorithm.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "Process.h"

namespace ProcessLib
{
struct ProcessData
{
    ProcessData(
        std::unique_ptr<NumLib::TimeStepAlgorithm>&& timestep_algorithm_,
        NumLib::NonlinearSolverTag const nonlinear_solver_tag_,
        NumLib::NonlinearSolverBase& nonlinear_solver_,
        std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit_,
        std::unique_ptr<NumLib::TimeDiscretization>&& time_disc_,
        int const process_id_, std::string&& process_name_, Process& process_)
        : timestep_algorithm(std::move(timestep_algorithm_)),
          timestep_previous(timestep_algorithm->begin()),
          timestep_current(timestep_previous),
          nonlinear_solver_tag(nonlinear_solver_tag_),
          nonlinear_solver(nonlinear_solver_),
          nonlinear_solver_status{true, 0},
          conv_crit(std::move(conv_crit_)),
          time_disc(std::move(time_disc_)),
          process_id(process_id_),
          process_name(std::move(process_name_)),
          process(process_)
    {
    }

    ProcessData(ProcessData&& pd) = delete;
    ProcessData& operator=(ProcessData const& pd) = delete;
    ProcessData& operator=(ProcessData&& pd) = delete;

    std::unique_ptr<NumLib::TimeStepAlgorithm> timestep_algorithm;
    NumLib::TimeStep timestep_previous;
    NumLib::TimeStep timestep_current;

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
    std::string const process_name;

    Process& process;
};

void setEquationSystem(ProcessData const& process_data);

}  // namespace ProcessLib
