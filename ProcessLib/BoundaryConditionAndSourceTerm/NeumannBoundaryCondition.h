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
#include "NeumannBoundaryConditionLocalAssembler.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
using NeumannBoundaryCondition =
    GenericNaturalBoundaryCondition<NeumannBoundaryConditionData,
                                    NeumannBoundaryConditionLocalAssembler>;

std::pair<std::string, std::optional<std::string>>
parseNeumannBoundaryCondition(BaseLib::ConfigTree const& config);

std::unique_ptr<NeumannBoundaryCondition> createNeumannBoundaryCondition(
    std::string const& parameter_name,
    std::optional<std::string> const& area_parameter_name,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
