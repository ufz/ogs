/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"

#include "ProcessLib/LIE/Common/FractureProperty.h"


namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace LIE
{

std::unique_ptr<BoundaryCondition>
createNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, bool is_axially_symmetric,
    unsigned const integration_order,
    unsigned const shapefunction_order,
    unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    FractureProperty const& fracture_prop);

}  // LIE
}  // ProcessLib
