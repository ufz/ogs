// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <string>
#include <vector>

namespace MeshLib
{

class Mesh;
class Node;
class Element;
}  // namespace MeshLib

namespace MeshToolsLib
{

/// Adds a layer to the mesh. If on_top is true, the layer is added on top,
/// if it is false, the layer is added at the bottom.
MeshLib::Mesh* addLayerToMesh(MeshLib::Mesh const& mesh, double const thickness,
                              std::string const& name, bool const on_top,
                              bool const copy_material_ids,
                              std::optional<int> const layer_id);

}  // namespace MeshToolsLib
