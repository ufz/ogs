/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>
#include <vector>

#include "MeshLib/Mesh.h"

namespace GeoLib
{
class GEOObjects;
}

namespace MeshGeoToolsLib
{
/// For each named geometry in the give geo_objects (defined on the given \c
/// mesh) constructs a mesh corresponding to the geometry with mappings to the
/// bulk mesh elements and nodes.
std::vector<std::unique_ptr<MeshLib::Mesh>>
constructAdditionalMeshesFromGeoObjects(GeoLib::GEOObjects const& geo_objects,
                                        MeshLib::Mesh const& mesh);

}  // namespace MeshGeoToolsLib
