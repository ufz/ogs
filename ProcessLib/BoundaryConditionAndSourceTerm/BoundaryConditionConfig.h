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

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
struct BoundaryConditionConfig final
{
    BaseLib::ConfigTree config;
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> boundary_meshes;
    std::optional<int> const component_id;
    bool compensate_non_equilibrium_initial_residuum = false;
};

}  // namespace ProcessLib
