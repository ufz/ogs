// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SourceTermCollection.h"

#include <range/v3/view/filter.hpp>

namespace ProcessLib
{
void SourceTermCollection::addSourceTermsForProcessVariables(
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    unsigned const integration_order,
    const MeshLib::Mesh& bulk_mesh)
{
    for (int variable_id = 0;
         variable_id < static_cast<int>(process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = process_variables[variable_id];
        auto sts =
            pv.createSourceTerms(dof_table, variable_id, integration_order,
                                 _parameters, process_variables, bulk_mesh);

        std::move(sts.begin(), sts.end(), std::back_inserter(_source_terms));
    }
}

void SourceTermCollection::integrate(const double t, GlobalVector const& x,
                                     GlobalVector& b, GlobalMatrix* jac) const
{
    // For parallel computing with DDC, a partition may not have source term
    // but a nullptr is assigned to its element in _source_terms.
    auto non_nullptr = [](std::unique_ptr<SourceTermBase> const& st)
    { return st != nullptr; };

    for (auto const& st : _source_terms | ranges::views::filter(non_nullptr))
    {
        st->integrate(t, x, b, jac);
    }
}

}  // namespace ProcessLib
