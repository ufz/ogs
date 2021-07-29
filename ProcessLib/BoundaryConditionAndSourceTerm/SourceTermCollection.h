/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/ProcessVariable.h"
#include "SourceTerm.h"

namespace ProcessLib
{
class SourceTermCollection final
{
public:
    explicit SourceTermCollection(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters)
        : _parameters(parameters)
    {
    }

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const;

    void addSourceTermsForProcessVariables(
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            process_variables,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        unsigned const integration_order);

private:
    std::vector<std::unique_ptr<SourceTerm>> _source_terms;
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        _parameters;
};

}  // namespace ProcessLib
