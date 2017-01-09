/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/IndexValueVector.h"
#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib
{
class BoundaryConditionCollection final
{
public:
    BoundaryConditionCollection(
        std::vector<std::unique_ptr<ParameterBase>> const& parameters)
        : _parameters(parameters)
    {
    }

    void applyNaturalBC(const double t, GlobalVector const& x, GlobalMatrix& K,
                        GlobalVector& b);

    std::vector<NumLib::IndexValueVector<GlobalIndexType>> const*
    getKnownSolutions(double const t) const
    {
        for (std::size_t i=0; i<_boundary_conditions.size(); ++i) {
            auto const& bc = *_boundary_conditions[i];
            auto& dirichlet_storage = _dirichlet_bcs[i];
            bc.getEssentialBCValues(t, dirichlet_storage);
        }
        return &_dirichlet_bcs;
    }

    void addBCsForProcessVariables(
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            process_variables,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        unsigned const integration_order);

    void preTimestep(const double t);

private:
    mutable std::vector<NumLib::IndexValueVector<GlobalIndexType>> _dirichlet_bcs;
    std::vector<std::unique_ptr<BoundaryCondition>> _boundary_conditions;
    std::vector<std::unique_ptr<ParameterBase>> const& _parameters;
};


}  // ProcessLib
