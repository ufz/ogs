/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 23, 2023, 2:22 PM
 */
#include "StaggeredCoupling.h"

#include <range/v3/view/filter.hpp>

#include "BaseLib/Error.h"
#include "BaseLib/RunTime.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/ODESolver/NonlinearSolver.h"

namespace NumLib
{
StaggeredCoupling::~StaggeredCoupling()
{
    for (auto* x : solutions_of_last_cpl_iteration_)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }
}

void StaggeredCoupling::initializeCoupledSolutions(
    std::vector<GlobalVector*> const& process_solutions)
{
    for (auto const* const x : process_solutions)
    {
        // Create a vector to store the solution of the last coupling iteration
        auto& x0 = NumLib::GlobalVectorProvider::provider.getVector(*x);
        MathLib::LinAlg::copy(*x, x0);

        // append a solution vector of suitable size
        solutions_of_last_cpl_iteration_.emplace_back(&x0);
    }
}

void StaggeredCoupling::setFirstIterationIndicator(
    std::vector<CouplingNodeVariant> const& coupling_nodes)
{
    auto is_regular_node = [](const CouplingNodeVariant& node)
    { return std::holds_alternative<CouplingNode>(node); };

    for (auto const& coupling_node :
         coupling_nodes | ranges::views::filter(is_regular_node))
    {
        std::get<CouplingNode>(coupling_node)
            .convergence_criterion->preFirstIteration();
    }
}

void StaggeredCoupling::resetCouplingConvergenceCriteria(
    std::vector<CouplingNodeVariant> const& coupling_nodes)
{
    auto is_regular_node = [](const CouplingNodeVariant& node)
    { return std::holds_alternative<CouplingNode>(node); };

    for (auto& coupling_node :
         coupling_nodes | ranges::views::filter(is_regular_node))
    {
        std::get<CouplingNode>(coupling_node).convergence_criterion->reset();
    }
}

bool StaggeredCoupling::checkCouplingConvergence(
    const bool convergence_of_last_process,
    CouplingNode const& coupling_node,
    GlobalVector const& x) const
{
    bool is_coupling_iteration_converged = convergence_of_last_process;

    auto& x_old = *solutions_of_last_cpl_iteration_[coupling_node.process_id];
    // Since x_old can be immediately refreshed after computing dx,
    // it is assigned with dx to save memory usage
    MathLib::LinAlg::axpy(x_old, -1.0, x);  // save dx = x - x_old to x_old.
    INFO(
        "------- Checking convergence criterion for coupled "
        "solution of process {:s} with ID {:d} -------",
        coupling_node.process_name, coupling_node.process_id);
    // Note: x_old stores dx
    coupling_node.convergence_criterion->checkDeltaX(x_old, x);

    is_coupling_iteration_converged =
        is_coupling_iteration_converged &&
        coupling_node.convergence_criterion->isSatisfied();

    return is_coupling_iteration_converged;
}

void StaggeredCoupling::updatePreviousSolution(int const process_id,
                                               GlobalVector const& x)
{
    auto& x_old = *solutions_of_last_cpl_iteration_[process_id];
    MathLib::LinAlg::copy(x, x_old);
}

}  // namespace NumLib
