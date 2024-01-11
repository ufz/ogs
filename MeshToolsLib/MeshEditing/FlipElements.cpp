/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FlipElements.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"

namespace MeshToolsLib
{
std::unique_ptr<MeshLib::Element> createFlippedElement(
    MeshLib::Element const& elem, std::vector<MeshLib::Node*> const& nodes)
{
    if (elem.getDimension() > 2)
    {
        return nullptr;
    }

    unsigned const n_nodes(elem.getNumberOfNodes());
    auto elem_nodes = std::make_unique<MeshLib::Node*[]>(n_nodes);
    for (unsigned i = 0; i < n_nodes; ++i)
    {
        elem_nodes[i] = nodes[elem.getNode(i)->getID()];
    }
    std::swap(elem_nodes[0], elem_nodes[1]);

    if (elem.getGeomType() == MeshLib::MeshElemType::LINE)
    {
        return std::make_unique<MeshLib::Line>(elem_nodes.release(),
                                               elem.getID());
    }
    if (elem.getGeomType() == MeshLib::MeshElemType::TRIANGLE)
    {
        return std::make_unique<MeshLib::Tri>(elem_nodes.release(),
                                              elem.getID());
    }
    if (elem.getGeomType() == MeshLib::MeshElemType::QUAD)
    {
        std::swap(elem_nodes[2], elem_nodes[3]);
        return std::make_unique<MeshLib::Quad>(elem_nodes.release(),
                                               elem.getID());
    }
    return nullptr;
}

std::unique_ptr<MeshLib::Mesh> createFlippedMesh(MeshLib::Mesh const& mesh)
{
    if (mesh.getDimension() > 2)
    {
        return nullptr;
    }

    std::vector<MeshLib::Node*> new_nodes(copyNodeVector(mesh.getNodes()));
    std::vector<MeshLib::Element*> const& elems(mesh.getElements());
    std::vector<MeshLib::Element*> new_elems;
    std::size_t n_elems(mesh.getNumberOfElements());
    new_elems.reserve(n_elems);

    for (std::size_t i = 0; i < n_elems; ++i)
    {
        new_elems.push_back(
            createFlippedElement(*elems[i], new_nodes).release());
    }

    return std::make_unique<MeshLib::Mesh>(
        "FlippedElementMesh", new_nodes, new_elems,
        true /* compute_element_neighbors */, mesh.getProperties());
}

}  // namespace MeshToolsLib
