/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 6, 2023, 4:21 PM
 */

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
