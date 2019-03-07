/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GenericNaturalBoundaryCondition.h"
#include "HCOpenBoundaryConditionLocalAssembler.h"
#include "MeshLib/PropertyVector.h"

namespace ProcessLib
{
using HCOpenBoundaryCondition =
    GenericNaturalBoundaryCondition<HCOpenBoundaryConditionData,
                                    HCOpenBoundaryConditionLocalAssembler>;

std::unique_ptr<HCOpenBoundaryCondition> createHCOpenBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const global_dim, Process const& process,
    unsigned const shapefunction_order);

}  // namespace ProcessLib
