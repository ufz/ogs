// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
        unsigned const integration_order,
        const MeshLib::Mesh& bulk_mesh);

private:
    std::vector<std::unique_ptr<SourceTermBase>> _source_terms;
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        _parameters;
};

}  // namespace ProcessLib
