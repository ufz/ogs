/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "QuadraticMeshGenerator.h"

#include <set>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"

/// Given an (linear) element divide all its edges by inserting a point in the
/// middle and return a new element.
template <typename QuadraticElement>
std::unique_ptr<QuadraticElement> convertLinearToQuadratic(
    MeshLib::Element const& e)
{
    int const n_all_nodes = QuadraticElement::n_all_nodes;
    int const n_base_nodes = QuadraticElement::n_base_nodes;
    assert(n_base_nodes == e.getNumberOfBaseNodes());

    // Copy base nodes of element to the quadratic element new nodes'.
    std::array<MeshLib::Node*, n_all_nodes> nodes{};
    for (int i = 0; i < n_base_nodes; i++)
    {
        nodes[i] = const_cast<MeshLib::Node*>(e.getNode(i));
    }

    // For each edge create a middle node.
    int const number_of_edges = e.getNumberOfEdges();
    for (int i = 0; i < number_of_edges; i++)
    {
        auto const& a = *e.getEdgeNode(i, 0);
        auto const& b = *e.getEdgeNode(i, 1);

        nodes[n_base_nodes + i] = new MeshLib::Node(
            (a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2);
    }

    return std::make_unique<QuadraticElement>(nodes, e.getID());
}

/// Special case for Quad-9 adding a centre node too.
template <>
std::unique_ptr<MeshLib::Quad9> convertLinearToQuadratic<MeshLib::Quad9>(
    MeshLib::Element const& e)
{
    int const n_all_nodes = MeshLib::Quad9::n_all_nodes;
    int const n_base_nodes = MeshLib::Quad9::n_base_nodes;
    assert(n_base_nodes == e.getNumberOfBaseNodes());

    // Copy base nodes of element to the quadratic element new nodes'.
    std::array<MeshLib::Node*, n_all_nodes> nodes{};
    for (int i = 0; i < n_base_nodes; i++)
    {
        nodes[i] = const_cast<MeshLib::Node*>(e.getNode(i));
    }

    // For each edge create a middle node.
    int const number_of_edges = e.getNumberOfEdges();
    for (int i = 0; i < number_of_edges; i++)
    {
        auto const& a = *e.getEdgeNode(i, 0);
        auto const& b = *e.getEdgeNode(i, 1);

        nodes[n_base_nodes + i] = new MeshLib::Node(
            (a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2);
    }

    // Compute the centre point coordinates.
    auto* centre_node = new MeshLib::Node(0, 0, 0);
    for (int i = 0; i < n_base_nodes; i++)
    {
        for (int d = 0; d < 3; d++)
        {
            (*centre_node)[d] += (*nodes[i])[d] / n_base_nodes;
        }
    }
    nodes[n_all_nodes - 1] = centre_node;

    return std::make_unique<MeshLib::Quad9>(nodes, e.getID());
}

/// Special case for Prism15
template <>
std::unique_ptr<MeshLib::Prism15> convertLinearToQuadratic<MeshLib::Prism15>(
    MeshLib::Element const& e)
{
    int const n_all_nodes = MeshLib::Prism15::n_all_nodes;
    int const n_base_nodes = MeshLib::Prism15::n_base_nodes;
    assert(n_base_nodes == e.getNumberOfBaseNodes());

    // Copy base nodes of element to the quadratic element new nodes'.
    std::array<MeshLib::Node*, n_all_nodes> nodes{};
    for (int i = 0; i < n_base_nodes; i++)
    {
        nodes[i] = const_cast<MeshLib::Node*>(e.getNode(i));
    }

    // For each edge create a middle node.
    int const number_of_edges = e.getNumberOfEdges();
    for (int i = 0; i < 3; i++)
    {
        auto const& a = *e.getEdgeNode(i, 0);
        auto const& b = *e.getEdgeNode(i, 1);

        nodes[n_base_nodes + i] = new MeshLib::Node(
            (a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2);
    }
    for (int i = 3; i < 6; i++)
    {
        auto const& a = *e.getEdgeNode(i + 3, 0);
        auto const& b = *e.getEdgeNode(i + 3, 1);

        nodes[n_base_nodes + i] = new MeshLib::Node(
            (a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2);
    }
    for (int i = 6; i < number_of_edges; i++)
    {
        auto const& a = *e.getEdgeNode(i - 3, 0);
        auto const& b = *e.getEdgeNode(i - 3, 1);

        nodes[n_base_nodes + i] = new MeshLib::Node(
            (a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2);
    }

    return std::make_unique<MeshLib::Prism15>(nodes, e.getID());
}

/// Return a new quadratic element corresponding to the linear element's type.
std::unique_ptr<MeshLib::Element> createQuadraticElement(
    MeshLib::Element const& e, bool const add_centre_node)
{
    if (e.getCellType() == MeshLib::CellType::LINE2)
    {
        return convertLinearToQuadratic<MeshLib::Line3>(e);
    }
    if (e.getCellType() == MeshLib::CellType::TRI3)
    {
        return convertLinearToQuadratic<MeshLib::Tri6>(e);
    }
    if (e.getCellType() == MeshLib::CellType::TET4)
    {
        return convertLinearToQuadratic<MeshLib::Tet10>(e);
    }
    if (e.getCellType() == MeshLib::CellType::QUAD4)
    {
        if (add_centre_node)
        {
            return convertLinearToQuadratic<MeshLib::Quad9>(e);
        }
        return convertLinearToQuadratic<MeshLib::Quad8>(e);
    }
    if (e.getCellType() == MeshLib::CellType::HEX8)
    {
        return convertLinearToQuadratic<MeshLib::Hex20>(e);
    }
    if (e.getCellType() == MeshLib::CellType::PRISM6)
    {
        return convertLinearToQuadratic<MeshLib::Prism15>(e);
    }
    if (e.getCellType() == MeshLib::CellType::PYRAMID5)
    {
        return convertLinearToQuadratic<MeshLib::Pyramid13>(e);
    }

    OGS_FATAL("Mesh element type {:s} is not supported",
              MeshLib::CellType2String(e.getCellType()));
}

struct nodeByCoordinatesComparator
{
    bool operator()(MeshLib::Node const* a, MeshLib::Node const* b) const
    {
        return *a < *b;
    }
};

namespace MeshToolsLib
{
std::unique_ptr<MeshLib::Mesh> createQuadraticOrderMesh(
    MeshLib::Mesh const& linear_mesh, bool const add_centre_node)
{
    // Clone the linear mesh nodes.
    auto quadratic_mesh_nodes = MeshLib::copyNodeVector(linear_mesh.getNodes());

    // Temporary container for unique quadratic nodes with O(log(n)) search.
    std::set<MeshLib::Node*, nodeByCoordinatesComparator> unique_nodes;

    // Create new elements with the quadratic nodes
    std::vector<MeshLib::Element*> quadratic_elements;
    auto const& linear_mesh_elements = linear_mesh.getElements();
    for (MeshLib::Element const* e : linear_mesh_elements)
    {
        auto quadratic_element = createQuadraticElement(*e, add_centre_node);

        // Replace the base nodes with cloned linear nodes.
        int const number_base_nodes = quadratic_element->getNumberOfBaseNodes();
        for (int i = 0; i < number_base_nodes; ++i)
        {
            quadratic_element->setNode(
                i, quadratic_mesh_nodes[getNodeIndex(*quadratic_element, i)]);
        }

        // Make the new (middle-edge) nodes unique.
        int const number_all_nodes = quadratic_element->getNumberOfNodes();
        for (int i = number_base_nodes; i < number_all_nodes; ++i)
        {
            MeshLib::Node* original_node =
                const_cast<MeshLib::Node*>(quadratic_element->getNode(i));

            auto it = unique_nodes.insert(original_node);
            if (!it.second)  // same node was already inserted before, no
                             // insertion
            {
                // Replace the element's node with the unique node.
                quadratic_element->setNode(i, *it.first);
                // And delete the original node
                delete original_node;
            }
        }

        quadratic_elements.push_back(quadratic_element.release());
    }

    // Add the unique quadratic nodes to the cloned linear nodes.
    quadratic_mesh_nodes.reserve(linear_mesh.getNodes().size() +
                                 unique_nodes.size());
    std::copy(unique_nodes.begin(), unique_nodes.end(),
              std::back_inserter(quadratic_mesh_nodes));

    return std::make_unique<MeshLib::Mesh>(
        linear_mesh.getName(), quadratic_mesh_nodes, quadratic_elements,
        true /* compute_element_neighbors */,
        linear_mesh.getProperties().excludeCopyProperties(
            std::vector<MeshLib::MeshItemType>(1,
                                               MeshLib::MeshItemType::Node)));
}

}  // namespace MeshToolsLib
