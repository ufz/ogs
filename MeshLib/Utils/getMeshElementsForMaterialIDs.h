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

#include <vector>

namespace MeshLib
{
class Element;
class Mesh;

std::vector<MeshLib::Element*> getMeshElementsForMaterialIDs(
    MeshLib::Mesh const& mesh, std::vector<int> const& selected_material_ids);
}  // namespace MeshLib
