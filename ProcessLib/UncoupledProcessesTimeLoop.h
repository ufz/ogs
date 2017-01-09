/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_UNCOUPLED_PROCESSES_TIMELOOP
#define PROCESSLIB_UNCOUPLED_PROCESSES_TIMELOOP

#include <memory>
#include <logog/include/logog.hpp>

#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/TimeStepping/Algorithms/ITimeStepAlgorithm.h"

#include "Output.h"
#include "Process.h"

namespace ProcessLib
{
struct SingleProcessData;

//! Time loop capable of time-integrating several uncoupled processes at once.
class UncoupledProcessesTimeLoop
{
public:
    explicit UncoupledProcessesTimeLoop(
        std::unique_ptr<NumLib::ITimeStepAlgorithm>&& timestepper,
        std::unique_ptr<Output>&& output,
        std::vector<std::unique_ptr<SingleProcessData>>&& per_process_data);

    bool loop();

    ~UncoupledProcessesTimeLoop();

private:
    std::vector<GlobalVector*> _process_solutions;
    std::unique_ptr<NumLib::ITimeStepAlgorithm> _timestepper;
    std::unique_ptr<Output> _output;
    std::vector<std::unique_ptr<SingleProcessData>> _per_process_data;
};

//! Builds an UncoupledProcessesTimeLoop from the given configuration.
std::unique_ptr<UncoupledProcessesTimeLoop> createUncoupledProcessesTimeLoop(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::map<std::string, std::unique_ptr<Process>> const&
        processes,
    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>> const&
        nonlinear_solvers);

}  // namespace ProcessLib

#endif  // PROCESSLIB_UNCOUPLED_PROCESSES_TIMELOOP
