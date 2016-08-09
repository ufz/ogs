/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_UNIFORMROBINBOUNDARYCONDITION_H
#define PROCESSLIB_UNIFORMROBINBOUNDARYCONDITION_H

#include "GenericNaturalBoundaryCondition.h"
#include "UniformRobinBoundaryConditionLocalAssembler.h"

namespace MeshGeoToolsLib
{
class BoundaryElementsSearcher;
}

namespace ProcessLib
{
using UniformRobinBoundaryCondition = GenericNaturalBoundaryCondition<
    UniformRobinBoundaryConditionData,
    UniformRobinBoundaryConditionLocalAssembler>;

/*! Creates a new uniform Robin boundary condition from the given data.
 *
 * The Robin boundary condition is given in the form
 * \f$ \partial f/\partial n = \alpha \cdot [ f(x) - f_0 ] \f$,
 * where the coefficients \f$ \alpha \f$ and \f$ f_0 \f$ are obtained from the
 * \c config.
 */
std::unique_ptr<UniformRobinBoundaryCondition>
createUniformRobinBoundaryCondition(
    BaseLib::ConfigTree const& config,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const global_dim);

}  // ProcessLib

#endif  // PROCESSLIB_UNIFORMROBINBOUNDARYCONDITION_H
