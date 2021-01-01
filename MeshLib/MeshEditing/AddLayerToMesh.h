/**
 * \file
 * \author Karsten Rink
 * \date   2016-01-18
 * \brief  Definition of AddLayerToMesh class
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <string>

namespace MeshLib
{

class Mesh;
class Node;
class Element;

/// Adds a layer to the mesh. If on_top is true, the layer is added on top,
/// if it is false, the layer is added at the bottom.
MeshLib::Mesh* addLayerToMesh(MeshLib::Mesh const& mesh, double thickness,
                              std::string const& name, bool on_top,
                              bool copy_material_ids);

} // end namespace MeshLib
