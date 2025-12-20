// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
