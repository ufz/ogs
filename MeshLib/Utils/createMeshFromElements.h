/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CREATEMESHFROMELEMENTS_H
#define CREATEMESHFROMELEMENTS_H

#include <vector>
#include <string>

#include "MeshLib/MeshEditing/MeshRevision.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
MeshLib::Mesh* createMeshFromElements(
    std::vector<MeshLib::Element*> const& elements,
    std::string const& new_mesh_name)
{
    std::vector<MeshLib::Element*> new_elements;
    new_elements.reserve(elements.size());
    for (auto const* e : elements)
    {
        new_elements.emplace_back(e->clone());
    }

   std::vector<MeshLib::Node*> new_nodes;
    for (auto * e : new_elements)
    {
        for (std::size_t k(0); k < e->getNumberOfNodes(); ++k)
        {
            new_nodes.push_back(new MeshLib::Node(*(e->getNodes()[k])));
            e->setNode(k, new_nodes.back());
        }
    }

    std::unique_ptr<MeshLib::Mesh> m(
        new MeshLib::Mesh("tmp", new_nodes, new_elements));
    MeshLib::MeshRevision mr(*m);
    return mr.collapseNodes(new_mesh_name,
                            std::numeric_limits<double>::epsilon());
}
}  // end namespace MeshLib

#endif
