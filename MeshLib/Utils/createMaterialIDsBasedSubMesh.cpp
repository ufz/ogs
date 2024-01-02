/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief Implementation of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
