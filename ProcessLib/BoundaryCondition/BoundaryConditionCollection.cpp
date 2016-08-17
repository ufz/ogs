/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundaryConditionCollection.h"

void initializeDirichletBCs(
    std::vector<std::unique_ptr<ProcessLib::BoundaryCondition>> const&
        boundary_conditions,
    std::vector<NumLib::IndexValueVector<GlobalIndexType>>&
        dirichlet_bcs)
{
    for (auto const& bc : boundary_conditions) {
        if (auto* dirichlet_bc =
                dynamic_cast<ProcessLib::DirichletBoundaryCondition*>(
                    bc.get())) {
            dirichlet_bcs.emplace_back(dirichlet_bc->getBCValues());
        }
    }
}

namespace ProcessLib
{
void BoundaryConditionCollection::apply(const double t, GlobalVector const& x,
                                        GlobalMatrix& K,
                                        GlobalVector& b)
{
    for (auto const& bc : _boundary_conditions)
        bc->apply(t, x, K, b);
}

void BoundaryConditionCollection::addBCsForProcessVariables(
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    unsigned const integration_order)
{
    for (int variable_id = 0;
         variable_id < static_cast<int>(process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = process_variables[variable_id];
        auto bcs = pv.createBoundaryConditions(dof_table, variable_id,
                                               integration_order, _parameters);

        std::move(bcs.begin(), bcs.end(),
                  std::back_inserter(_boundary_conditions));
    }

    initializeDirichletBCs(_boundary_conditions, _dirichlet_bcs);
}
}
