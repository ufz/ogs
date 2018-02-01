/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/ProcessVariable.h"
#include "ProcessLib/SourceTerms/NodalSourceTerm.h"

namespace ProcessLib
{
class SourceTermCollection final
{
public:
    SourceTermCollection(
        std::vector<std::unique_ptr<ParameterBase>> const& parameters)
        : _parameters(parameters)
    {
    }

    void integrateNodalSourceTerms(const double t, GlobalVector& b) const;

    void addSourceTermsForProcessVariables(
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            process_variables,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        unsigned const integration_order);

private:
    std::vector<std::unique_ptr<NodalSourceTerm>> _source_terms;
    std::vector<std::unique_ptr<ParameterBase>> const& _parameters;
};

}  // ProcessLib
