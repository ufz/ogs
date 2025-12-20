// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "GenericNaturalBoundaryCondition.h"
#include "HCNonAdvectiveFreeComponentFlowBoundaryConditionLocalAssembler.h"
#include "MeshLib/PropertyVector.h"

namespace ProcessLib
{
using HCNonAdvectiveFreeComponentFlowBoundaryCondition =
    GenericNaturalBoundaryCondition<
        HCNonAdvectiveFreeComponentFlowBoundaryConditionData,
        HCNonAdvectiveFreeComponentFlowBoundaryConditionLocalAssembler>;

std::string parseHCNonAdvectiveFreeComponentFlowBoundaryCondition(
    BaseLib::ConfigTree const& config);

std::unique_ptr<HCNonAdvectiveFreeComponentFlowBoundaryCondition>
createHCNonAdvectiveFreeComponentFlowBoundaryCondition(
    std::string const& boundary_permeability_name, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const global_dim, Process const& process,
    unsigned const shapefunction_order);

}  // namespace ProcessLib
