/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
