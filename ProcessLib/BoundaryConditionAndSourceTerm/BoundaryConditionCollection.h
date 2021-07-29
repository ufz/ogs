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

#include "BoundaryCondition.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib
{
class BoundaryConditionCollection final
{
public:
    explicit BoundaryConditionCollection(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters)
        : _parameters(parameters)
    {
    }

    void applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                        int const process_id, GlobalMatrix& K, GlobalVector& b,
                        GlobalMatrix* Jac) const;

    std::vector<NumLib::IndexValueVector<GlobalIndexType>> const*
    getKnownSolutions(double const t, GlobalVector const& x) const
    {
        auto const n_bcs = _boundary_conditions.size();
        for (std::size_t i = 0; i < n_bcs; ++i)
        {
            auto const& bc = *_boundary_conditions[i];
            auto& dirichlet_storage = _dirichlet_bcs[i];
            bc.getEssentialBCValues(t, x, dirichlet_storage);
        }
        return &_dirichlet_bcs;
    }

    void addBCsForProcessVariables(
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            process_variables,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        unsigned const integration_order, Process const& process);

    void addBoundaryCondition(std::unique_ptr<BoundaryCondition>&& bc)
    {
        _boundary_conditions.push_back(std::move(bc));
    }

    void preTimestep(const double t, std::vector<GlobalVector*> const& x,
                     int const process_id) const
    {
        for (auto const& bc_ptr : _boundary_conditions)
        {
            bc_ptr->preTimestep(t, x, process_id);
        }
    }

    void postTimestep(const double t, std::vector<GlobalVector*> const& x,
                      int const process_id) const
    {
        for (auto const& bc_ptr : _boundary_conditions)
        {
            bc_ptr->postTimestep(t, x, process_id);
        }
    }

private:
    mutable std::vector<NumLib::IndexValueVector<GlobalIndexType>>
        _dirichlet_bcs;
    std::vector<std::unique_ptr<BoundaryCondition>> _boundary_conditions;
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        _parameters;
};

}  // namespace ProcessLib
