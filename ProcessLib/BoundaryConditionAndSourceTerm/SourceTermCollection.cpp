/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SourceTermCollection.h"

#include <range/v3/view/filter.hpp>

namespace ProcessLib
{
void SourceTermCollection::addSourceTermsForProcessVariables(
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
        auto sts =
            pv.createSourceTerms(dof_table, variable_id, integration_order,
                                 _parameters, process_variables);

        std::move(sts.begin(), sts.end(), std::back_inserter(_source_terms));
    }
}

void SourceTermCollection::integrate(const double t, GlobalVector const& x,
                                     GlobalVector& b, GlobalMatrix* jac) const
{
    // For parallel computing with DDC, a partition may not have source term
    // but a nullptr is assigned to its element in _source_terms.
    auto non_nullptr = [](std::unique_ptr<SourceTerm> const& st)
    { return st != nullptr; };

    for (auto const& st : _source_terms | ranges::views::filter(non_nullptr))
    {
        st->integrate(t, x, b, jac);
    }
}

}  // namespace ProcessLib
