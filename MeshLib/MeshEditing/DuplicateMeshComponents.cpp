/**
 * \file
 * \author Karsten Rink
 * \date   2014-03-25
 * \brief  Implementation of Duplicate functions.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DuplicateMeshComponents.h"

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
std::vector<Node*> copyNodeVector(const std::vector<Node*>& nodes)
{
    const std::size_t nNodes(nodes.size());
    std::vector<Node*> new_nodes;
    new_nodes.reserve(nNodes);
    for (std::size_t k = 0; k < nNodes; ++k)
    {
        new_nodes.push_back(new Node(nodes[k]->getCoords(), new_nodes.size()));
    }
    return new_nodes;
}

std::vector<Element*> copyElementVector(
    std::vector<Element*> const& elements,
    std::vector<Node*> const& new_nodes,
    std::vector<std::size_t> const* const node_id_map)
{
    std::vector<Element*> new_elements;
    new_elements.reserve(elements.size());
    std::transform(elements.begin(), elements.end(),
                   std::back_inserter(new_elements),
                   [&new_nodes, &node_id_map](auto const& element) {
                       return copyElement(element, new_nodes, node_id_map);
                   });
    return new_elements;
}

/// Copies an element without change, using the nodes vector from the result
/// mesh.
template <typename E>
Element* copyElement(Element const* const element,
                     const std::vector<Node*>& nodes,
                     std::vector<std::size_t> const* const id_map)
{
    unsigned const number_of_element_nodes(element->getNumberOfNodes());
    auto** new_nodes = new Node*[number_of_element_nodes];
    if (id_map)
    {
        for (unsigned i = 0; i < number_of_element_nodes; ++i)
        {
            new_nodes[i] = nodes[(*id_map)[element->getNode(i)->getID()]];
        }
    }
    else
    {
        for (unsigned i = 0; i < number_of_element_nodes; ++i)
        {
            new_nodes[i] = nodes[element->getNode(i)->getID()];
        }
    }
    return new E(new_nodes);
}

Element* copyElement(Element const* const element,
                     const std::vector<Node*>& nodes,
                     std::vector<std::size_t> const* const id_map)
{
    switch (element->getCellType())
    {
        case CellType::LINE2:
            return copyElement<Line>(element, nodes, id_map);
        case CellType::LINE3:
            return copyElement<Line3>(element, nodes, id_map);
        case CellType::TRI3:
            return copyElement<Tri>(element, nodes, id_map);
        case CellType::TRI6:
            return copyElement<Tri6>(element, nodes, id_map);
        case CellType::QUAD4:
            return copyElement<Quad>(element, nodes, id_map);
        case CellType::QUAD8:
            return copyElement<Quad8>(element, nodes, id_map);
        case CellType::QUAD9:
            return copyElement<Quad9>(element, nodes, id_map);
        case CellType::TET4:
            return copyElement<Tet>(element, nodes, id_map);
        case CellType::TET10:
            return copyElement<Tet10>(element, nodes, id_map);
        case CellType::HEX8:
            return copyElement<Hex>(element, nodes, id_map);
        case CellType::HEX20:
            return copyElement<Hex20>(element, nodes, id_map);
        case CellType::PYRAMID5:
            return copyElement<Pyramid>(element, nodes, id_map);
        case CellType::PYRAMID13:
            return copyElement<Pyramid13>(element, nodes, id_map);
        case CellType::PRISM6:
            return copyElement<Prism>(element, nodes, id_map);
        case CellType::PRISM15:
            return copyElement<Prism15>(element, nodes, id_map);
        default:
        {
            ERR("Error: Unknown cell type.");
            return nullptr;
        }
    }
}

std::vector<Element*> cloneElements(std::vector<Element*> const& elements)
{
    std::vector<Element*> cloned_elements;
    cloned_elements.reserve(elements.size());
    std::transform(begin(elements), end(elements),
                   std::back_inserter(cloned_elements),
                   [](Element* const e) { return e->clone(); });
    return cloned_elements;
}

}  // namespace MeshLib
