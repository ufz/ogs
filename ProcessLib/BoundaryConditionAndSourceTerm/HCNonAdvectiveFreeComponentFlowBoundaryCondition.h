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
#include "HCNonAdvectiveFreeComponentFlowBoundaryConditionLocalAssembler.h"
#include "MeshLib/PropertyVector.h"

namespace ProcessLib
{
using HCNonAdvectiveFreeComponentFlowBoundaryCondition =
    GenericNaturalBoundaryCondition<
        HCNonAdvectiveFreeComponentFlowBoundaryConditionData,
        HCNonAdvectiveFreeComponentFlowBoundaryConditionLocalAssembler>;

std::unique_ptr<HCNonAdvectiveFreeComponentFlowBoundaryCondition>
createHCNonAdvectiveFreeComponentFlowBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const global_dim, Process const& process,
    unsigned const shapefunction_order);

}  // namespace ProcessLib
