// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

VariableDependentNeumannConfig parseVariableDependentNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config);

std::unique_ptr<VariableDependentNeumannBoundaryCondition>
createVariableDependentNeumannBoundaryCondition(
    VariableDependentNeumannConfig const& coefficients,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
