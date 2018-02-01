/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SourceTermCollection.h"

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
        auto sts = pv.createSourceTerms(dof_table, variable_id,
                                        integration_order, _parameters);

        std::move(sts.begin(), sts.end(), std::back_inserter(_source_terms));
    }
}

void SourceTermCollection::integrateNodalSourceTerms(const double t,
                                                     GlobalVector& b) const
{
    for (auto const& st : _source_terms)
        st->integrateNodalSourceTerm(t, b);
}

}
