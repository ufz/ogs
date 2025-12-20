// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

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
