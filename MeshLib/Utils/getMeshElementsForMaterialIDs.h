// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace MeshLib
{
class Element;
class Mesh;

std::vector<MeshLib::Element*> getMeshElementsForMaterialIDs(
    MeshLib::Mesh const& mesh, std::vector<int> const& selected_material_ids);
}  // namespace MeshLib
