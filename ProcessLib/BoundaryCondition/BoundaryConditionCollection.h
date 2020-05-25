/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/IndexValueVector.h"
#include "ProcessLib/ProcessVariable.h"

#include "BoundaryCondition.h"

namespace ProcessLib
{
class BoundaryConditionCollection final
{
public:
    explicit BoundaryConditionCollection(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters)
        : parameters_(parameters)
    {
    }

    void applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                        int const process_id, GlobalMatrix& K, GlobalVector& b,
                        GlobalMatrix* Jac);

    std::vector<NumLib::IndexValueVector<GlobalIndexType>> const*
    getKnownSolutions(double const t, GlobalVector const& x) const
    {
        auto const n_bcs = boundary_conditions_.size();
        for (std::size_t i=0; i<n_bcs; ++i) {
            auto const& bc = *boundary_conditions_[i];
            auto& dirichlet_storage = dirichlet_bcs_[i];
            bc.getEssentialBCValues(t, x, dirichlet_storage);
        }
        return &dirichlet_bcs_;
    }

    void addBCsForProcessVariables(
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            process_variables,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        unsigned const integration_order, Process const& process);

    void addBoundaryCondition(std::unique_ptr<BoundaryCondition>&& bc)
    {
        boundary_conditions_.push_back(std::move(bc));
    }

    void preTimestep(const double t, std::vector<GlobalVector*> const& x,
                     int const process_id)
    {
        for (auto const& bc_ptr : boundary_conditions_)
        {
            bc_ptr->preTimestep(t, x, process_id);
        }
    }

private:
    mutable std::vector<NumLib::IndexValueVector<GlobalIndexType>> dirichlet_bcs_;
    std::vector<std::unique_ptr<BoundaryCondition>> boundary_conditions_;
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters_;
};

}  // namespace ProcessLib
