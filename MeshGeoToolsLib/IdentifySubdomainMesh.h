/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

namespace MeshGeoToolsLib
{
class SearchLength;
class MeshNodeSearcher;
}
namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{
/// Geometrically finds nodes and elements of the subdomain mesh in the bulk
/// mesh, and updates or verifies the corresponding bulk_node_ids and
/// bulk_element_ids properties.
///
/// In case of unexpected results OGS_FATAL is called.
void identifySubdomainMesh(MeshLib::Mesh& subdomain_mesh,
                           MeshLib::Mesh const& bulk_mesh,
                           MeshNodeSearcher const& mesh_node_searcher,
                           bool const force_overwrite = false);
}  // namespace MeshGeoToolsLib
