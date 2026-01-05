// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "GenericNaturalBoundaryCondition.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MeshLib/PropertyVector.h"
#include "WellboreCompensateNeumannBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
using WellboreCompensateNeumannBoundaryCondition =
    GenericNaturalBoundaryCondition<
        WellboreCompensateNeumannBoundaryConditionData,
        WellboreCompensateNeumannBoundaryConditionLocalAssembler>;

WellboreCompensateCoefficients parseWellboreCompensateNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config);

std::unique_ptr<WellboreCompensateNeumannBoundaryCondition>
createWellboreCompensateNeumannBoundaryCondition(
    WellboreCompensateCoefficients const& coefficients,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace ProcessLib
