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
#include <unordered_map>
#include <typeinfo>
#include <typeindex>

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

/// Time loop capable of time-integrating several processes at once.
/// TODO: Rename to, e.g., TimeLoop, since it is not for purely uncoupled stuff
/// anymore.
class UncoupledProcessesTimeLoop
{
public:
    explicit UncoupledProcessesTimeLoop(
        std::unique_ptr<Output>&& output,
        std::vector<std::unique_ptr<SingleProcessData>>&& per_process_data,
        const unsigned global_coupling_max_iterations,
        std::unique_ptr<NumLib::ConvergenceCriterion>&&
            global_coupling_conv_crit,
        const double start_time, const double end_time);

    bool loop();

    ~UncoupledProcessesTimeLoop();

    /**
     *  This function fills the vector of solutions of coupled processes of
     *  processes, _solutions_of_coupled_processes, and initializes the vector
     * of
     *  solutions of the previous coupling iteration,
     *  _solutions_of_last_cpl_iteration.
     *
     *  \return a boolean value as a flag to indicate there should be a coupling
     *          among processes or not.
     */
    bool setCoupledSolutions();

private:
    std::vector<GlobalVector*> _process_solutions;
    std::unique_ptr<Output> _output;
    std::vector<std::unique_ptr<SingleProcessData>> _per_process_data;

    const double _start_time;
    const double _end_time;

    /// Maximum iterations of the global coupling.
    const unsigned _global_coupling_max_iterations;
    /// Convergence criteria of the global coupling iterations.
    std::unique_ptr<NumLib::ConvergenceCriterion> _global_coupling_conv_crit;

    /**
     *  Vector of solutions of coupled processes of processes.
     *  Each vector element stores the references of the solution vectors
     *  (stored in _process_solutions) of the coupled processes of a process.
     */
    std::vector<std::unordered_map<std::type_index, GlobalVector const&>>
        _solutions_of_coupled_processes;

    /// Solutions of the previous coupling iteration for the convergence
    /// criteria of the coupling iteration.
    std::vector<GlobalVector*> _solutions_of_last_cpl_iteration;

    /**
     * \brief Member to solver non coupled systems of equations, which can be
     *        a single system of equations, or several systems of equations
     *        without any dependency among the different systems.
     *
     * @param t           Current time
     * @param dt          Time step size
     * @param timestep_id Index of the time step
     * @return            true:  if all nonlinear solvers convergence.
     *                    false: if any of nonlinear solvers divergences.
     */
    bool solveUncoupledEquationSystems(const double t, const double dt,
                                       const std::size_t timestep_id);

    /**
     * \brief Member to solver coupled systems of equations by the staggered
     *        scheme.
     *
     * @param t           Current time
     * @param dt          Time step size
     * @param timestep_id Index of the time step
     * @return            true:   if all nonlinear solvers convergence.
     *                    false:  if any of nonlinear solvers divergences.
     */
    bool solveCoupledEquationSystemsByStaggeredScheme(
        const double t, const double dt, const std::size_t timestep_id);
};

//! Builds an UncoupledProcessesTimeLoop from the given configuration.
std::unique_ptr<UncoupledProcessesTimeLoop> createUncoupledProcessesTimeLoop(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::map<std::string, std::unique_ptr<Process>> const& processes,
    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>> const&
        nonlinear_solvers);

}  // namespace ProcessLib
