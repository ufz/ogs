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

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
struct BoundaryConditionConfig final
{
    BoundaryConditionConfig(BaseLib::ConfigTree&& config_,
                            MeshLib::Mesh const& mesh_,
                            std::optional<int> const component_id_)
        : config(std::move(config_)),
          boundary_mesh(mesh_),
          component_id(component_id_)
    {
    }

    BoundaryConditionConfig(BoundaryConditionConfig&& other) = default;

    BaseLib::ConfigTree config;
    MeshLib::Mesh const& boundary_mesh;
    std::optional<int> const component_id;
};

}  // namespace ProcessLib
