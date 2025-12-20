// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
