/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Parameter/Parameter.h"
#include "GenericNaturalBoundaryCondition.h"
#include "NeumannBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
using NeumannBoundaryCondition = GenericNaturalBoundaryCondition<
    Parameter<double> const&, NeumannBoundaryConditionLocalAssembler>;

std::unique_ptr<NeumannBoundaryCondition> createNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters);

}  // ProcessLib
