/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_NEUMANNBOUNDARYCONDITION_H
#define PROCESSLIB_NEUMANNBOUNDARYCONDITION_H

#include "ProcessLib/Parameter/Parameter.h"
#include "GenericNaturalBoundaryCondition.h"
#include "NeumannBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
using NeumannBoundaryCondition = GenericNaturalBoundaryCondition<
    Parameter<double> const&, NeumannBoundaryConditionLocalAssembler>;

std::unique_ptr<NeumannBoundaryCondition> createNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, bool is_axially_symmetric,
    unsigned const integration_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters);

}  // ProcessLib

#endif  // PROCESSLIB_NEUMANNBOUNDARYCONDITION_H
