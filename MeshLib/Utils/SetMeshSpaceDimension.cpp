/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 5, 2021, 12:46 PM
 */

#include "SetMeshSpaceDimension.h"

#include <range/v3/algorithm/for_each.hpp>

#include "GetSpaceDimension.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
void setMeshSpaceDimension(std::vector<std::unique_ptr<Mesh>> const& meshes)
{
    // Get the space dimension from the bulk mesh:
    auto const d = getSpaceDimension(meshes[0]->getNodes());
    for (auto const& mesh : meshes)
    {
        ranges::for_each(mesh->getElements(),
                         [d](Element* const e) { e->space_dimension_ = d; });
    }
}
};  // namespace MeshLib
