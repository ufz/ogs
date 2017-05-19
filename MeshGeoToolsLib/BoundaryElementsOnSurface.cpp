/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "BoundaryElementsOnSurface.h"

#include "GeoLib/Surface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshSearch/ElementSearch.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{
BoundaryElementsOnSurface::BoundaryElementsOnSurface(
    MeshLib::Mesh const& mesh, MeshNodeSearcher const& mshNodeSearcher,
    GeoLib::Surface const& sfc)
    : _mesh(mesh), _sfc(sfc)
{
    // search elements near the surface
    auto node_ids_on_sfc = mshNodeSearcher.getMeshNodeIDsAlongSurface(sfc);
    MeshLib::ElementSearch es(_mesh);
    es.searchByNodeIDs(node_ids_on_sfc);
    auto &ele_ids_near_sfc = es.getSearchedElementIDs();

    // get a list of faces made of the nodes
    for (auto ele_id : ele_ids_near_sfc) {
        auto* e = _mesh.getElement(ele_id);
        // skip internal elements
        if (!e->isBoundaryElement())
            continue;
        // find faces on surface
        for (unsigned i=0; i<e->getNumberOfFaces(); i++) {
            auto* face = e->getFace(i);
            // check
            std::size_t cnt_match = 0;
            for (std::size_t j=0; j<face->getNumberOfBaseNodes(); j++) {
                if (std::find(node_ids_on_sfc.begin(), node_ids_on_sfc.end(), face->getNodeIndex(j)) != node_ids_on_sfc.end())
                    cnt_match++;
                else
                    break;
            }
            // update the list
            if (cnt_match==face->getNumberOfBaseNodes())
                _boundary_elements.push_back(const_cast<MeshLib::Element*>(face));
            else
                delete face;
        }
    }
}

BoundaryElementsOnSurface::~BoundaryElementsOnSurface()
{
    for (auto p : _boundary_elements)
        delete p;
}

} // end namespace MeshGeoToolsLib

