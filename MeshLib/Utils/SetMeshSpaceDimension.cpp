/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 5, 2021, 12:46 PM
 */

#include "SetMeshSpaceDimension.h"

#include "GetSpaceDimension.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
void setMeshSpaceDimension(std::vector<std::unique_ptr<Mesh>> const& meshes)
{
    // Get the space dimension from the bulk mesh:
    auto const space_dimension = getSpaceDimension(meshes[0]->getNodes());
    for (auto& mesh : meshes)
    {
        auto elements = mesh->getElements();
        for (auto element : elements)
        {
            element->space_dimension_ = space_dimension;
        }
    }
}
};  // namespace MeshLib
