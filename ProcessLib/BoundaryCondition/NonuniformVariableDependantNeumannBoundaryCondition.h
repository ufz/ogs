/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GenericNonuniformNaturalBoundaryCondition.h"
#include "MeshLib/PropertyVector.h"
#include "NonuniformVariableDependantNeumannBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
using NonuniformVariableDependantNeumannBoundaryCondition =
    GenericNonuniformNaturalBoundaryCondition<
        NonuniformVariableDependantNeumannBoundaryConditionData,
        NonuniformVariableDependantNeumannBoundaryConditionLocalAssembler>;

std::unique_ptr<NonuniformVariableDependantNeumannBoundaryCondition>
createNonuniformVariableDependantNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, const MeshLib::Mesh& bulk_mesh);

}  // namespace ProcessLib
