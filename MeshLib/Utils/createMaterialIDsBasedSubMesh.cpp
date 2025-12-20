// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshLib/Utils/createMaterialIDsBasedSubMesh.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"
#include "MeshLib/Utils/createMeshFromElementSelection.h"
#include "MeshLib/Utils/getMeshElementsForMaterialIDs.h"

namespace MeshLib
{
std::unique_ptr<MeshLib::Mesh> createMaterialIDsBasedSubMesh(
    MeshLib::Mesh const& mesh, std::vector<int> const& material_ids,
    std::string const& name_for_created_mesh)
{
    auto const elements =
        MeshLib::getMeshElementsForMaterialIDs(mesh, material_ids);
    return MeshLib::createMeshFromElementSelection(
        name_for_created_mesh, MeshLib::cloneElements(elements));
}
}  // namespace MeshLib
