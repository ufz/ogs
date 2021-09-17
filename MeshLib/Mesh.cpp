/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief Implementation of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Mesh.h"

#include <memory>
#include <unordered_map>
#include <utility>

#include "BaseLib/RunTime.h"
#include "Elements/Element.h"
#include "Elements/Hex.h"
#include "Elements/Prism.h"
#include "Elements/Pyramid.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Tri.h"

/// Mesh counter used to uniquely identify meshes by id.
static std::size_t global_mesh_counter = 0;

namespace MeshLib
{
Mesh::Mesh(std::string name,
           std::vector<Node*>
               nodes,
           std::vector<Element*>
               elements,
           Properties const& properties)
    : _id(global_mesh_counter++),
      _mesh_dimension(0),
      _edge_length(std::numeric_limits<double>::max(), 0),
      _node_distance(std::numeric_limits<double>::max(), 0),
      _name(std::move(name)),
      _nodes(std::move(nodes)),
      _elements(std::move(elements)),
      _properties(properties)
{
    this->resetNodeIDs();
    this->resetElementIDs();
    this->setDimension();
    this->setElementsConnectedToNodes();
    this->setElementNeighbors();

    this->calcEdgeLengthRange();
}

Mesh::Mesh(const Mesh& mesh)
    : _id(global_mesh_counter++),
      _mesh_dimension(mesh.getDimension()),
      _edge_length(mesh._edge_length.first, mesh._edge_length.second),
      _node_distance(mesh._node_distance.first, mesh._node_distance.second),
      _name(mesh.getName()),
      _nodes(mesh.getNumberOfNodes()),
      _elements(mesh.getNumberOfElements()),
      _properties(mesh._properties)
{
    const std::vector<Node*>& nodes(mesh.getNodes());
    const std::size_t nNodes(nodes.size());
    for (unsigned i = 0; i < nNodes; ++i)
    {
        _nodes[i] = new Node(*nodes[i]);
    }

    const std::vector<Element*>& elements(mesh.getElements());
    const std::size_t nElements(elements.size());
    for (unsigned i = 0; i < nElements; ++i)
    {
        const std::size_t nElemNodes = elements[i]->getNumberOfNodes();
        _elements[i] = elements[i]->clone();
        for (unsigned j = 0; j < nElemNodes; ++j)
        {
            _elements[i]->setNode(j, _nodes[elements[i]->getNode(j)->getID()]);
        }
    }

    if (_mesh_dimension == 0)
    {
        this->setDimension();
    }
    this->setElementsConnectedToNodes();
    this->setElementNeighbors();
}

Mesh::~Mesh()
{
    const std::size_t nElements(_elements.size());
    for (std::size_t i = 0; i < nElements; ++i)
    {
        delete _elements[i];
    }

    const std::size_t nNodes(_nodes.size());
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        delete _nodes[i];
    }
}

void Mesh::addElement(Element* elem)
{
    _elements.push_back(elem);

    // add element information to nodes
    unsigned nNodes(elem->getNumberOfNodes());
    for (unsigned i = 0; i < nNodes; ++i)
    {
        elem->getNode(i)->addElement(elem);
    }
}

void Mesh::resetNodeIDs()
{
    const std::size_t nNodes(_nodes.size());
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        _nodes[i]->setID(i);
    }
}

void Mesh::resetElementIDs()
{
    const std::size_t nElements(this->_elements.size());
    for (unsigned i = 0; i < nElements; ++i)
    {
        _elements[i]->setID(i);
    }
}

void Mesh::setDimension()
{
    const std::size_t nElements(_elements.size());
    for (unsigned i = 0; i < nElements; ++i)
    {
        if (_elements[i]->getDimension() > _mesh_dimension)
        {
            _mesh_dimension = _elements[i]->getDimension();
        }
    }
}

void Mesh::setElementsConnectedToNodes()
{
    for (auto& element : _elements)
    {
        const unsigned nNodes(element->getNumberOfNodes());
        for (unsigned j = 0; j < nNodes; ++j)
        {
            element->getNode(j)->addElement(element);
        }
    }
}

void Mesh::calcEdgeLengthRange()
{
    const std::size_t nElems(getNumberOfElements());
    for (std::size_t i = 0; i < nElems; ++i)
    {
        auto const& [min_length, max_length] =
            computeSqrEdgeLengthRange(*_elements[i]);
        _edge_length.first = std::min(_edge_length.first, min_length);
        _edge_length.second = std::max(_edge_length.second, max_length);
    }
    _edge_length.first = std::sqrt(_edge_length.first);
    _edge_length.second = std::sqrt(_edge_length.second);
}

void Mesh::setElementNeighbors()
{
    std::vector<Element*> neighbors;
    for (auto element : _elements)
    {
        // create vector with all elements connected to current element
        // (includes lots of doubles!)
        const std::size_t nNodes(element->getNumberOfBaseNodes());
        for (unsigned n(0); n < nNodes; ++n)
        {
            std::vector<Element*> const& conn_elems(
                (element->getNode(n)->getElements()));
            neighbors.insert(neighbors.end(), conn_elems.begin(),
                             conn_elems.end());
        }
        std::sort(neighbors.begin(), neighbors.end());
        auto const neighbors_new_end =
            std::unique(neighbors.begin(), neighbors.end());

        for (auto neighbor = neighbors.begin(); neighbor != neighbors_new_end;
             ++neighbor)
        {
            std::optional<unsigned> const opposite_face_id =
                element->addNeighbor(*neighbor);
            if (opposite_face_id)
            {
                (*neighbor)->setNeighbor(element, *opposite_face_id);
            }
        }
        neighbors.clear();
    }
}

std::size_t Mesh::getNumberOfBaseNodes() const
{
    return std::count_if(begin(_nodes), end(_nodes),
                         [](auto const* const node)
                         { return isBaseNode(*node); });
}

bool Mesh::hasNonlinearElement() const
{
    return std::any_of(
        std::begin(_elements), std::end(_elements), [](Element const* const e) {
            return e->getNumberOfNodes() != e->getNumberOfBaseNodes();
        });
}

void scaleMeshPropertyVector(MeshLib::Mesh& mesh,
                             std::string const& property_name,
                             double factor)
{
    if (!mesh.getProperties().existsPropertyVector<double>(property_name))
    {
        WARN("Did not find PropertyVector '{:s}' for scaling.", property_name);
        return;
    }
    auto& pv = *mesh.getProperties().getPropertyVector<double>(property_name);
    std::transform(pv.begin(), pv.end(), pv.begin(),
                   [factor](auto const& v) { return v * factor; });
}

PropertyVector<int> const* materialIDs(Mesh const& mesh)
{
    auto const& properties = mesh.getProperties();
    return properties.existsPropertyVector<int>("MaterialIDs",
                                                MeshLib::MeshItemType::Cell, 1)
               ? properties.getPropertyVector<int>(
                     "MaterialIDs", MeshLib::MeshItemType::Cell, 1)
               : nullptr;
}

std::unique_ptr<MeshLib::Mesh> createMeshFromElementSelection(
    std::string mesh_name, std::vector<MeshLib::Element*> const& elements)
{
    DBUG("Found {:d} elements in the mesh", elements.size());

    // Store bulk element ids for each of the new elements.
    std::vector<std::size_t> bulk_element_ids;
    bulk_element_ids.reserve(elements.size());
    std::transform(begin(elements), end(elements),
                   std::back_inserter(bulk_element_ids),
                   [&](auto const& e) { return e->getID(); });

    // original node pointers to newly created nodes.
    std::unordered_map<const MeshLib::Node*, MeshLib::Node*> nodes_map;
    nodes_map.reserve(
        elements.size());  // There will be at least one node per element.

    for (auto& e : elements)
    {
        // For each node find a cloned node in map or create if there is none.
        unsigned const n_nodes = e->getNumberOfNodes();
        for (unsigned i = 0; i < n_nodes; ++i)
        {
            const MeshLib::Node* n = e->getNode(i);
            auto const it = nodes_map.find(n);
            if (it == nodes_map.end())
            {
                auto new_node_in_map = nodes_map[n] = new MeshLib::Node(*n);
                e->setNode(i, new_node_in_map);
            }
            else
            {
                e->setNode(i, it->second);
            }
        }
    }

    // Copy the unique nodes pointers.
    std::vector<MeshLib::Node*> element_nodes;
    element_nodes.reserve(nodes_map.size());
    std::transform(begin(nodes_map), end(nodes_map),
                   std::back_inserter(element_nodes),
                   [](auto const& pair) { return pair.second; });

    // Store bulk node ids for each of the new nodes.
    std::vector<std::size_t> bulk_node_ids;
    bulk_node_ids.reserve(nodes_map.size());
    std::transform(begin(nodes_map), end(nodes_map),
                   std::back_inserter(bulk_node_ids),
                   [](auto const& pair) { return pair.first->getID(); });

    auto mesh = std::make_unique<MeshLib::Mesh>(
        std::move(mesh_name), std::move(element_nodes), std::move(elements));
    assert(mesh != nullptr);

    addPropertyToMesh(*mesh, "bulk_element_ids", MeshLib::MeshItemType::Cell, 1,
                      bulk_element_ids);
    addPropertyToMesh(*mesh, "bulk_node_ids", MeshLib::MeshItemType::Node, 1,
                      bulk_node_ids);

    return mesh;
}

std::vector<std::vector<Node*>> calculateNodesConnectedByElements(
    Mesh const& mesh)
{
    std::vector<std::vector<Node*>> nodes_connected_by_elements;
    auto const& nodes = mesh.getNodes();
    nodes_connected_by_elements.resize(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        auto& adjacent_nodes = nodes_connected_by_elements[i];
        auto const* node = nodes[i];

        // Get all elements, to which this node is connected.
        auto const& connected_elements = node->getElements();

        // And collect all elements' nodes.
        for (Element const* const element : connected_elements)
        {
            Node* const* const single_elem_nodes = element->getNodes();
            std::size_t const nnodes = element->getNumberOfNodes();
            for (std::size_t n = 0; n < nnodes; n++)
            {
                adjacent_nodes.push_back(single_elem_nodes[n]);
            }
        }

        // Make nodes unique and sorted by their ids.
        // This relies on the node's id being equivalent to it's address.
        std::sort(adjacent_nodes.begin(), adjacent_nodes.end(),
                  [](Node* a, Node* b) { return a->getID() < b->getID(); });
        auto const last =
            std::unique(adjacent_nodes.begin(), adjacent_nodes.end());
        adjacent_nodes.erase(last, adjacent_nodes.end());
    }
    return nodes_connected_by_elements;
}

}  // namespace MeshLib
