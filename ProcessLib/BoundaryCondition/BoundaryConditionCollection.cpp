/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundaryConditionCollection.h"

namespace ProcessLib
{
void BoundaryConditionCollection::applyNaturalBC(
    const double t, std::vector<GlobalVector*> const& x, int const process_id,
    GlobalMatrix& K, GlobalVector& b, GlobalMatrix* Jac)
{
    for (auto const& bc : boundary_conditions_)
    {
        bc->applyNaturalBC(t, x, process_id, K, b, Jac);
    }
}

void BoundaryConditionCollection::addBCsForProcessVariables(
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    unsigned const integration_order, Process const& process)
{
    for (int variable_id = 0;
         variable_id < static_cast<int>(process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = process_variables[variable_id];
        auto bcs = pv.createBoundaryConditions(
            dof_table, variable_id, integration_order, parameters_, process);

        std::move(bcs.begin(), bcs.end(),
                  std::back_inserter(boundary_conditions_));
    }

    // For each BC there will be storage for Dirichlet BC. This storage will be
    // uninitialized by default, and has to be filled by the respective BC
    // object if needed.
    dirichlet_bcs_.resize(boundary_conditions_.size());
}
}  // namespace ProcessLib
