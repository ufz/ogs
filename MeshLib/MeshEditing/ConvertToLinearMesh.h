/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>

#include "MeshLib/Mesh.h"

namespace MeshLib
{

/// Converts a non-linear mesh to a linear mesh. All the mesh properties will
/// be copied except for entries for non-linear nodes.
std::unique_ptr<MeshLib::Mesh> convertToLinearMesh(
    const MeshLib::Mesh& mesh, const std::string& new_mesh_name);

} // end namespace MeshLib
