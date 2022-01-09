/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FlipElements.h"

#include "DuplicateMeshComponents.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
std::unique_ptr<Element> createFlippedElement(Element const& elem,
                                              std::vector<Node*> const& nodes)
{
    if (elem.getDimension() > 2)
    {
        return nullptr;
    }

    unsigned const n_nodes(elem.getNumberOfNodes());
    auto elem_nodes = std::make_unique<Node*[]>(n_nodes);
    for (unsigned i = 0; i < n_nodes; ++i)
    {
        elem_nodes[i] = nodes[elem.getNode(i)->getID()];
    }
    std::swap(elem_nodes[0], elem_nodes[1]);

    if (elem.getGeomType() == MeshElemType::LINE)
    {
        return std::make_unique<Line>(elem_nodes.release(), elem.getID());
    }
    if (elem.getGeomType() == MeshElemType::TRIANGLE)
    {
        return std::make_unique<Tri>(elem_nodes.release(), elem.getID());
    }
    if (elem.getGeomType() == MeshElemType::QUAD)
    {
        std::swap(elem_nodes[2], elem_nodes[3]);
        return std::make_unique<Quad>(elem_nodes.release(), elem.getID());
    }
    return nullptr;
}

std::unique_ptr<Mesh> createFlippedMesh(Mesh const& mesh)
{
    if (mesh.getDimension() > 2)
    {
        return nullptr;
    }

    std::vector<Node*> new_nodes(copyNodeVector(mesh.getNodes()));
    std::vector<Element*> const& elems(mesh.getElements());
    std::vector<Element*> new_elems;
    std::size_t n_elems(mesh.getNumberOfElements());
    new_elems.reserve(n_elems);

    for (std::size_t i = 0; i < n_elems; ++i)
    {
        new_elems.push_back(
            createFlippedElement(*elems[i], new_nodes).release());
    }

    return std::make_unique<Mesh>("FlippedElementMesh", new_nodes, new_elems,
                                  mesh.getProperties());
}

}  // end namespace MeshLib
