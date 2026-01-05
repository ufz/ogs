// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace MeshLib
{
class Mesh;

}

namespace MeshToolsLib
{

void zeroMeshFieldDataByMaterialIDs(
    MeshLib::Mesh& mesh, std::vector<int> const& selected_material_ids);
}  // namespace MeshToolsLib
