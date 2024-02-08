/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 20, 2023, 5:09 PM
 */

#pragma once

#include <functional>
#include <memory>
#include <tuple>
#include <variant>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/ODESolver/NonlinearSolverStatus.h"

namespace NumLib
{
class ConvergenceCriterion;

/// Information of a coupling node
struct CouplingNode
{
    std::string process_name;
    std::unique_ptr<NumLib::ConvergenceCriterion> convergence_criterion;
    int max_iterations;
    int process_id;
};

struct RootCouplingNode
{
    int max_iterations;
    using CouplingNodeVariant = std::variant<CouplingNode, RootCouplingNode>;
    std::vector<CouplingNodeVariant> sub_coupling_nodes;
};

/// A class designed to manage the solution of coupled equations using the
/// staggered method.
class StaggeredCoupling
{
    template <typename ProcessData, typename Output>
    using ProcessSolver = std::function<NumLib::NonlinearSolverStatus(
        std::vector<GlobalVector*>& /*xs*/,
        std::vector<GlobalVector*> const& /*xs_prev*/,
        std::size_t const /*timestep*/, double const /*t*/,
        double const /*delta_t*/, ProcessData const& /*process_data*/,
        std::vector<Output> const& /*outputs*/)>;

    using CouplingNodeVariant = std::variant<CouplingNode, RootCouplingNode>;

public:
    StaggeredCoupling(const int global_coupling_max_iterations,
                      std::vector<CouplingNodeVariant>&& coupling_nodes)
        : global_coupling_max_iterations_(global_coupling_max_iterations),
          coupling_nodes_(std::move(coupling_nodes))
    {
    }

    ~StaggeredCoupling();

    /**
     * This function fills the vector of solutions of coupled processes of
     * processes, solutions_of_coupled_processes_, and initializes the
     * vector of solutions of the previous coupling iteration,
     * _solutions_of_last_cpl_iteration.
     */
    void initializeCoupledSolutions(
        std::vector<GlobalVector*> const& process_solutions);

    /**
     * It solves the equations of all coupled processes by the staggered method,
     * and it returns nonlinear solver status.
     */
    template <typename ProcessData, typename Output>
    NumLib::NonlinearSolverStatus execute(
        const double t, const double dt, const std::size_t timestep_id,
        std::vector<GlobalVector*>& process_solutions,
        std::vector<GlobalVector*> const& process_solutions_prev,
        std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
        std::vector<Output> const& outputs,
        ProcessSolver<ProcessData, Output> const&
            solve_one_time_step_one_process);

private:
    /// Maximum iteration number of the coupling loop of the staggered scheme.
    const int global_coupling_max_iterations_;

    /** Coupling graph for the staggered scheme. It looks like:
     *
     *    x ... x    o        o
     *              / \      / \
     *             /   \    /   \
     *            x ... x  x ... x
     *
     *   where x represents a coupling node, and o represents a dummy root
     *   coupling node.
     */
    std::vector<CouplingNodeVariant> coupling_nodes_;

    /**
     * It solves the equation of each coupling process (treated as a coupling
     * node) via function recursion, and it returns nonlinear solver status, and
     * an indicator of coupling iteration convergence.
     *
     */
    template <typename ProcessData, typename Output>
    std::tuple<NumLib::NonlinearSolverStatus, bool> executeConcrete(
        std::vector<CouplingNodeVariant>& coupling_nodes,
        const int max_iterations, const double t, const double dt,
        const std::size_t timestep_id,
        std::vector<GlobalVector*>& process_solutions,
        std::vector<GlobalVector*> const& process_solutions_prev,
        std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
        std::vector<Output> const& outputs,
        ProcessSolver<ProcessData, Output> const&
            solve_one_time_step_one_process);

    template <typename ProcessData, typename Output>
    std::tuple<NumLib::NonlinearSolverStatus, bool> executeSubCoupling(
        CouplingNodeVariant& coupling_node, const double t, const double dt,
        const std::size_t timestep_id,
        std::vector<GlobalVector*>& process_solutions,
        std::vector<GlobalVector*> const& process_solutions_prev,
        std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
        std::vector<Output> const& outputs,
        ProcessSolver<ProcessData, Output> const&
            solve_one_time_step_one_process);

    template <typename ProcessData, typename Output>
    NumLib::NonlinearSolverStatus executeSingleIteration(
        int const global_coupling_iteration,
        CouplingNode const& regular_coupling_node, const double t,
        const double dt, const std::size_t timestep_id,
        std::vector<GlobalVector*>& process_solutions,
        std::vector<GlobalVector*> const& process_solutions_prev,
        std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
        std::vector<Output> const& outputs,
        ProcessSolver<ProcessData, Output> const&
            solve_one_time_step_one_process);

    /// Solutions of the previous coupling iteration for the convergence
    /// criteria of the coupling iteration.
    std::vector<GlobalVector*> solutions_of_last_cpl_iteration_;

    /// Set the indicator of the first staggered coupling iteration be true.
    void setFirstIterationIndicator(
        std::vector<CouplingNodeVariant> const& coupling_nodes);

    void resetCouplingConvergenceCriteria(
        std::vector<CouplingNodeVariant> const& coupling_nodes);

    bool checkCouplingConvergence(const bool convergence_of_last_process,
                                  CouplingNode const& coupling_node,
                                  GlobalVector const& x) const;

    void updatePreviousSolution(int const process_id, GlobalVector const& x);
};

}  // namespace NumLib

#include "StaggeredCoupling-impl.h"
