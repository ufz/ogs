/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_UNIFORMNEUMANNBOUNDARYCONDITION_H
#define PROCESSLIB_UNIFORMNEUMANNBOUNDARYCONDITION_H

#include "GenericNaturalBoundaryCondition.h"
#include "UniformNeumannBoundaryConditionLocalAssembler.h"

namespace MeshGeoToolsLib
{
class BoundaryElementsSearcher;
}

namespace ProcessLib
{
using UniformNeumannBoundaryCondition = GenericNaturalBoundaryCondition<
    double, UniformNeumannBoundaryConditionLocalAssembler>;

std::unique_ptr<UniformNeumannBoundaryCondition>
createUniformNeumannBoundaryCondition(
    BoundaryConditionConfig const& config,
    MeshGeoToolsLib::BoundaryElementsSearcher& searcher,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    unsigned const integration_order, unsigned const global_dim);

}  // ProcessLib

#endif  // PROCESSLIB_UNIFORMNEUMANNBOUNDARYCONDITION_H
