/*
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshGeoToolsLib/convertMeshNodesToGeometry.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{
void convertMeshNodesToGeometry(MeshLib::Mesh const& mesh,
                                std::vector<std::size_t> const& node_ids,
                                std::string& geo_name,
                                GeoLib::GEOObjects& geometry_sets)
{
    std::vector<MeshLib::Node*> const& nodes(mesh.getNodes());
    auto pnts = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    std::map<std::string, std::size_t>* pnt_names(
        new std::map<std::string, std::size_t>);
    std::size_t cnt(0);
    for (std::size_t id : node_ids)
    {
        pnts->push_back(new GeoLib::Point(*(nodes[id]), cnt));
        pnt_names->insert(std::pair<std::string, std::size_t>(
            "node_id-"+std::to_string(nodes[id]->getID()), cnt));
        cnt++;
    }

    geometry_sets.addPointVec(std::move(pnts), geo_name, pnt_names);
}

}  // end namespace MeshGeoToolsLib
