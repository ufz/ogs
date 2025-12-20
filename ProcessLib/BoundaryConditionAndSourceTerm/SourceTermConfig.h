// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
struct SourceTermConfig final
{
    SourceTermConfig(BaseLib::ConfigTree&& config_,
                     MeshLib::Mesh const& mesh_,
                     int component_id_)
        : config(std::move(config_)), mesh(mesh_), component_id(component_id_)
    {
    }

    SourceTermConfig(SourceTermConfig&& other)
        : config(std::move(other.config)),
          mesh(other.mesh),
          component_id(other.component_id)
    {
    }

    BaseLib::ConfigTree config;
    MeshLib::Mesh const& mesh;
    int component_id;
};

}  // namespace ProcessLib
