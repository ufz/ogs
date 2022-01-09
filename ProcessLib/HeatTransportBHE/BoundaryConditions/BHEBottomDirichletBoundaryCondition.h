/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/IndexValueVector.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/BoundaryCondition.h"

namespace ProcessLib::HeatTransportBHE
{
class BHEBottomDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    explicit BHEBottomDirichletBoundaryCondition(
        std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices)
        : _in_out_global_indices(std::move(in_out_global_indices))
    {
    }

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    std::pair<GlobalIndexType, GlobalIndexType> const _in_out_global_indices;
};

std::unique_ptr<BHEBottomDirichletBoundaryCondition>
createBHEBottomDirichletBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices);
}  // namespace ProcessLib::HeatTransportBHE
