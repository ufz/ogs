/**
 * \file
 * \author Karsten Rink
 * \date   2016-01-18
 * \brief  Definition of AddLayerToMesh class
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
