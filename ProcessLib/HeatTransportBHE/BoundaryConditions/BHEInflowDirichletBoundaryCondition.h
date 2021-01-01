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

#include "NumLib/IndexValueVector.h"
#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"

namespace ProcessLib::HeatTransportBHE
{
template <typename BHEUpdateCallback>
class BHEInflowDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    BHEInflowDirichletBoundaryCondition(
        std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
        BHEUpdateCallback bhe_update_callback)
        : _in_out_global_indices(std::move(in_out_global_indices)),
          _bhe_update_callback(bhe_update_callback)
    {
    }

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override
    {
        bc_values.ids.resize(1);
        bc_values.values.resize(1);

        bc_values.ids[0] = _in_out_global_indices.first;
        // here call the corresponding BHE functions
        auto const T_out = x[_in_out_global_indices.second];
        bc_values.values[0] = _bhe_update_callback(T_out, t);
    }

private:
    std::pair<GlobalIndexType, GlobalIndexType> const _in_out_global_indices;
    BHEUpdateCallback _bhe_update_callback;
};

template <typename BHEUpdateCallback>
std::unique_ptr<BHEInflowDirichletBoundaryCondition<BHEUpdateCallback>>
createBHEInflowDirichletBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    BHEUpdateCallback bhe_update_callback)
{
    DBUG("Constructing BHEInflowDirichletBoundaryCondition.");

    // In case of partitioned mesh the boundary could be empty, i.e. there is no
    // boundary condition.
#ifdef USE_PETSC
    // For this special boundary condition the boundary condition is not empty
    // if the global indices are non-negative.
    if (in_out_global_indices.first < 0 && in_out_global_indices.second < 0)
    {
        return nullptr;
    }
    // If only one of the global indices (in or out) is negative the
    // implementation is not valid.
    if (in_out_global_indices.first < 0 || in_out_global_indices.second < 0)
    {
        OGS_FATAL(
            "The partition cuts the BHE into two independent parts. This "
            "behaviour is not implemented.");
    }
#endif  // USE_PETSC

    return std::make_unique<
        BHEInflowDirichletBoundaryCondition<BHEUpdateCallback>>(
        std::move(in_out_global_indices), bhe_update_callback);
}
}  // namespace ProcessLib::HeatTransportBHE
