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

#include "getMeshElementsForMaterialIDs.h"

#include <range/v3/algorithm/contains.hpp>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
std::vector<MeshLib::Element*> getMeshElementsForMaterialIDs(
    MeshLib::Mesh const& mesh, std::vector<int> const& selected_material_ids)
{
    auto const material_ids = *materialIDs(mesh);
    auto const& elements = mesh.getElements();
    std::vector<MeshLib::Element*> selected_elements;

    for (std::size_t i = 0; i < material_ids.size(); ++i)
    {
        if (ranges::contains(selected_material_ids, material_ids[i]))
        {
            selected_elements.push_back(elements[i]);
        }
    }
    return selected_elements;
}
}  // namespace MeshLib
