// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>

#include "MeshLib/Mesh.h"

namespace MeshToolsLib
{
/// Converts a non-linear mesh to a linear mesh. All the mesh properties will
/// be copied except for entries for non-linear nodes.
std::unique_ptr<MeshLib::Mesh> convertToLinearMesh(
    const MeshLib::Mesh& mesh, const std::string& new_mesh_name);

}  // namespace MeshToolsLib
