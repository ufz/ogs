// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
#include <vector>

namespace MeshLib
{
class Mesh;

std::unique_ptr<MeshLib::Mesh> createMaterialIDsBasedSubMesh(
    MeshLib::Mesh const& mesh, std::vector<int> const& material_ids,
    std::string const& name_for_created_mesh);

}  // namespace MeshLib
