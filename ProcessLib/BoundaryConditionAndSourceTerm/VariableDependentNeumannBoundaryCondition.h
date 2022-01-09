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

#include "GenericNaturalBoundaryCondition.h"
#include "MeshLib/PropertyVector.h"
#include "VariableDependentNeumannBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
using VariableDependentNeumannBoundaryCondition =
    GenericNaturalBoundaryCondition<
        VariableDependentNeumannBoundaryConditionData,
        VariableDependentNeumannBoundaryConditionLocalAssembler>;

std::unique_ptr<VariableDependentNeumannBoundaryCondition>
createVariableDependentNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
