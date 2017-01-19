/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <logog/include/logog.hpp>

#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/TimeStepping/Algorithms/ITimeStepAlgorithm.h"

#include "Output.h"
#include "Process.h"

namespace NumLib
{
class ConvergenceCriterion;
}

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
        std::vector<std::unique_ptr<SingleProcessData>>&& per_process_data,
        const unsigned global_coupling_max_iterations,
        std::unique_ptr<NumLib::ConvergenceCriterion>&& glb_coupling_conv_crit);

    bool loop();

    ~UncoupledProcessesTimeLoop();

    bool setCoupledSolutions();

private:
    std::vector<GlobalVector*> _process_solutions;
    std::unique_ptr<NumLib::ITimeStepAlgorithm> _timestepper;
    std::unique_ptr<Output> _output;
    std::vector<std::unique_ptr<SingleProcessData>> _per_process_data;

    /// Maximum iterations of the global coupling.
    const unsigned _global_coupling_max_iterations;
    /// Convergence criteria of the global coupling iterations.
    std::unique_ptr<NumLib::ConvergenceCriterion> _global_coupling_conv_crit;
    std::vector<std::map<ProcessType, GlobalVector const*>>
            _solutions_of_coupled_processes;
    /// Solutions of the previous coupling iteration.
    std::vector<GlobalVector*> _solutions_of_last_cpl_iteration;
};

//! Builds an UncoupledProcessesTimeLoop from the given configuration.
std::unique_ptr<UncoupledProcessesTimeLoop> createUncoupledProcessesTimeLoop(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::map<std::string, std::unique_ptr<Process>> const&
        processes,
    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>> const&
        nonlinear_solvers);

}  // namespace ProcessLib
