/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include "QuadraticMeshGenerator.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Node.h"

/// Given an (linear) element divide all its edges by inserting a point in the
/// middle and return a new element.
template <typename QuadraticElement>
std::unique_ptr<QuadraticElement> convertLinearToQuadratic(
    MeshLib::Element const& e)
{
    auto const n_all_nodes = QuadraticElement::n_all_nodes;
    auto const n_base_nodes = QuadraticElement::n_base_nodes;
    assert(n_base_nodes == e.getNumberOfBaseNodes());

    // Copy base nodes of element to the quadratic element new nodes'.
    std::array<MeshLib::Node*, n_all_nodes> nodes;
    for (int i = 0; i < n_base_nodes; i++)
        nodes[i] = const_cast<MeshLib::Node*>(e.getNode(i));

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

/// Return a new quadratic element corresponding to the linear element's type.
std::unique_ptr<MeshLib::Element> createQuadraticElement(
    MeshLib::Element const& e)
{
    if (e.getCellType() == MeshLib::CellType::LINE2)
    {
        return convertLinearToQuadratic<MeshLib::Line3>(e);
    }
    if (e.getCellType() == MeshLib::CellType::TRI3)
    {
        return convertLinearToQuadratic<MeshLib::Tri6>(e);
    }
    if (e.getCellType() == MeshLib::CellType::QUAD4)
    {
        return convertLinearToQuadratic<MeshLib::Quad8>(e);
    }
    if (e.getCellType() == MeshLib::CellType::HEX8)
    {
        return convertLinearToQuadratic<MeshLib::Hex20>(e);
    }

    OGS_FATAL("Mesh element type %s is not supported",
              MeshLib::CellType2String(e.getCellType()).c_str());
}

struct nodeByCoordinatesComparator
{
    bool operator()(MeshLib::Node* a, MeshLib::Node* b) const
    {
        return *a < *b;
    }
};

namespace MeshLib
{
std::unique_ptr<Mesh> createQuadraticOrderMesh(Mesh const& linear_mesh)
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
        auto quadratic_element = createQuadraticElement(*e);

        // Replace the base nodes with cloned linear nodes.
        int const number_base_nodes = quadratic_element->getNumberOfBaseNodes();
        for (int i = 0; i < number_base_nodes; ++i)
        {
            quadratic_element->setNode(
                i, quadratic_mesh_nodes[quadratic_element->getNodeIndex(i)]);
        }

        // Make the new (middle-edge) nodes unique.
        int const number_all_nodes = quadratic_element->getNumberOfNodes();
        for (int i = number_base_nodes; i < number_all_nodes; ++i)
        {
            Node* original_node =
                const_cast<Node*>(quadratic_element->getNode(i));

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
        linear_mesh.getProperties().excludeCopyProperties(
            std::vector<MeshLib::MeshItemType>(1, MeshLib::MeshItemType::Node)),
        linear_mesh.getNumberOfNodes());
}

}  // namespace MeshLib
