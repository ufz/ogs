/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLIB_CONVERRTTOLINEARMESH_H_
#define MESHLIB_CONVERRTTOLINEARMESH_H_

#include <memory>
#include <string>

#include "MeshLib/Mesh.h"

namespace MeshLib
{

/// Converts a non-linear mesh to a linear meash. All the mesh properties will
/// be copied except for entries for non-linear nodes.
std::unique_ptr<MeshLib::Mesh> convertToLinearMesh(
    const MeshLib::Mesh& mesh, const std::string& new_mesh_name);

} // end namespace MeshLib

#endif //MESHLIB_CONVERRTTOLINEARMESH_H_
