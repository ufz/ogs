/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEBottomDirichletBoundaryCondition.h"
#include "BaseLib/Error.h"

namespace ProcessLib::HeatTransportBHE
{
void BHEBottomDirichletBoundaryCondition::getEssentialBCValues(
    const double /*t*/, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    bc_values.ids.resize(1);
    bc_values.values.resize(1);

    bc_values.ids[0] = _in_out_global_indices.second;
    // here, the outflow temperature is always
    // the same as the inflow temperature
    // get the inflow temperature from here.
    bc_values.values[0] = x[_in_out_global_indices.first];
}

std::unique_ptr<BHEBottomDirichletBoundaryCondition>
createBHEBottomDirichletBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices)
{
    DBUG("Constructing BHEBottomDirichletBoundaryCondition.");

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

    return std::make_unique<BHEBottomDirichletBoundaryCondition>(
        std::move(in_out_global_indices));
}
}  // namespace ProcessLib::HeatTransportBHE
