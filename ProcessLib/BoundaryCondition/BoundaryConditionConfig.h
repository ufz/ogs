/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
                            boost::optional<int> const component_id_)
        : config(std::move(config_)),
          boundary_mesh(mesh_),
          component_id(component_id_)
    {
    }

    BoundaryConditionConfig(BoundaryConditionConfig&& other)
        : config(std::move(other.config)),
          boundary_mesh(other.boundary_mesh),
          component_id(other.component_id)
    {
    }

    BaseLib::ConfigTree config;
    MeshLib::Mesh const& boundary_mesh;
    boost::optional<int> const component_id;
};

}  // ProcessLib
