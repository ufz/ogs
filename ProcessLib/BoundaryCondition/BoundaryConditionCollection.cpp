/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundaryConditionCollection.h"

namespace ProcessLib
{
void BoundaryConditionCollection::applyNaturalBC(const double t,
                                                 GlobalVector const& x,
                                                 GlobalMatrix& K,
                                                 GlobalVector& b)
{
    for (auto const& bc : _boundary_conditions)
        bc->applyNaturalBC(t, x, K, b);
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
            dof_table, variable_id, integration_order, _parameters, process);

        std::move(bcs.begin(), bcs.end(),
                  std::back_inserter(_boundary_conditions));
    }

    // For each BC there will be storage for Dirichlet BC. This storage will be
    // uninitialized by default, and has to be filled by the respective BC
    // object if needed.
    _dirichlet_bcs.resize(_boundary_conditions.size());
}
}
