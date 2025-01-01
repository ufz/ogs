/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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

std::unique_ptr<WellboreCompensateNeumannBoundaryCondition>
createWellboreCompensateNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace ProcessLib
