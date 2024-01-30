/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief Implementation of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Mesh.h"

#include <memory>
#include <range/v3/algorithm/contains.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/map.hpp>
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
#include "Utils/DuplicateMeshComponents.h"
#include "Utils/getMeshElementsForMaterialIDs.h"

/// Mesh counter used to uniquely identify meshes by id.
static std::size_t global_mesh_counter = 0;

namespace MeshLib
{
using namespace ranges;

std::vector<std::vector<Element const*>> findElementsConnectedToNodes(
    Mesh const& mesh)
{
    std::vector<std::vector<Element const*>> elements_connected_to_nodes;
    auto const& nodes = mesh.getNodes();
    elements_connected_to_nodes.resize(nodes.size());

    for (auto const* element : mesh.getElements())
    {
        for (auto const node_id : element->nodes() | views::ids)
        {
            elements_connected_to_nodes[node_id].push_back(element);
        }
    }
    return elements_connected_to_nodes;
}

Mesh::Mesh(std::string name,
           std::vector<Node*>
               nodes,
           std::vector<Element*>
               elements,
           bool const compute_element_neighbors,
           Properties const& properties)
    : _id(global_mesh_counter++),
      _mesh_dimension(0),
      _node_distance(std::numeric_limits<double>::max(), 0),
      _name(std::move(name)),
      _nodes(std::move(nodes)),
      _elements(std::move(elements)),
      _properties(properties),
      _compute_element_neighbors(compute_element_neighbors)
{
    this->resetNodeIDs();
    this->resetElementIDs();
    this->setDimension();

    _elements_connected_to_nodes = findElementsConnectedToNodes(*this);

    if (_compute_element_neighbors)
    {
        this->setElementNeighbors();
    }
}

Mesh::Mesh(const Mesh& mesh)
    : _id(global_mesh_counter++),
      _mesh_dimension(mesh.getDimension()),
      _node_distance(mesh._node_distance.first, mesh._node_distance.second),
      _name(mesh.getName()),
      _nodes(mesh.getNumberOfNodes()),
      _elements(mesh.getNumberOfElements()),
      _properties(mesh._properties),
      _compute_element_neighbors(mesh._compute_element_neighbors)
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
        _elements[i] = elements[i]->clone();
        for (auto const& [j, node_id] :
             elements[i]->nodes() | views::ids | ranges::views::enumerate)
        {
            _elements[i]->setNode(static_cast<unsigned>(j), _nodes[node_id]);
        }
    }

    if (_mesh_dimension == 0)
    {
        this->setDimension();
    }

    _elements_connected_to_nodes = findElementsConnectedToNodes(*this);
    if (_compute_element_neighbors)
    {
        this->setElementNeighbors();
    }
}

Mesh::Mesh(Mesh&& mesh) = default;

void Mesh::shallowClean()
{
    _elements.clear();
    _nodes.clear();
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

std::pair<double, double> minMaxEdgeLength(
    std::vector<Element*> const& elements)
{
    auto min_max = [](auto const a, auto const b) -> std::pair<double, double> {
        return {std::min(a.first, b.first), std::max(a.second, b.second)};
    };

    using limits = std::numeric_limits<double>;
    auto const bounds = ranges::accumulate(
        elements, std::pair{limits::infinity(), -limits::infinity()}, min_max,
        [](Element const* const e) { return computeSqrEdgeLengthRange(*e); });

    return {std::sqrt(bounds.first), std::sqrt(bounds.second)};
}

void Mesh::setElementNeighbors()
{
    std::vector<Element const*> neighbors;
    for (auto element : _elements)
    {
        // create vector with all elements connected to current element
        // (includes lots of doubles!)
        const std::size_t nNodes(element->getNumberOfBaseNodes());
        for (unsigned n(0); n < nNodes; ++n)
        {
            auto const& conn_elems(
                _elements_connected_to_nodes[element->getNode(n)->getID()]);
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
                element->addNeighbor(const_cast<Element*>(*neighbor));
            if (opposite_face_id)
            {
                const_cast<Element*>(*neighbor)->setNeighbor(element,
                                                             *opposite_face_id);
            }
        }
        neighbors.clear();
    }
}

std::size_t Mesh::computeNumberOfBaseNodes() const
{
    return std::count_if(begin(_nodes), end(_nodes),
                         [this](auto const* const node) {
                             return isBaseNode(
                                 *node,
                                 _elements_connected_to_nodes[node->getID()]);
                         });
}

bool Mesh::hasNonlinearElement() const
{
    return std::any_of(
        std::begin(_elements), std::end(_elements),
        [](Element const* const e)
        { return e->getNumberOfNodes() != e->getNumberOfBaseNodes(); });
}

std::vector<MeshLib::Element const*> const& Mesh::getElementsConnectedToNode(
    std::size_t const node_id) const
{
    return _elements_connected_to_nodes[node_id];
}

std::vector<MeshLib::Element const*> const& Mesh::getElementsConnectedToNode(
    Node const& node) const
{
    return _elements_connected_to_nodes[node.getID()];
}

PropertyVector<int> const* materialIDs(Mesh const& mesh)
{
    auto const& properties = mesh.getProperties();
    if (properties.existsPropertyVector<int>("MaterialIDs",
                                             MeshLib::MeshItemType::Cell, 1))
    {
        return properties.getPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    }
    if (properties.hasPropertyVector("MaterialIDs"))
    {
        WARN(
            "The 'MaterialIDs' mesh property exists but is either of wrong "
            "type (must be int), or it is not defined on element / cell data.");
    }
    return nullptr;
}

PropertyVector<int>* materialIDs(Mesh& mesh)
{
    return const_cast<PropertyVector<int>*>(
        MeshLib::materialIDs(std::as_const(mesh)));
}

PropertyVector<std::size_t> const* bulkNodeIDs(Mesh const& mesh)
{
    auto const& properties = mesh.getProperties();
    return properties.getPropertyVector<std::size_t>(
        MeshLib::getBulkIDString(MeshLib::MeshItemType::Node),
        MeshLib::MeshItemType::Node, 1);
}

PropertyVector<std::size_t> const* bulkElementIDs(Mesh const& mesh)
{
    auto const& properties = mesh.getProperties();
    return properties.getPropertyVector<std::size_t>(
        MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell),
        MeshLib::MeshItemType::Cell, 1);
}

std::vector<std::vector<Node*>> calculateNodesConnectedByElements(
    Mesh const& mesh)
{
    auto const elements_connected_to_nodes = findElementsConnectedToNodes(mesh);

    std::vector<std::vector<Node*>> nodes_connected_by_elements;
    auto const& nodes = mesh.getNodes();
    nodes_connected_by_elements.resize(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        auto& adjacent_nodes = nodes_connected_by_elements[i];
        auto const node_id = nodes[i]->getID();

        // Get all elements, to which this node is connected.
        auto const& connected_elements = elements_connected_to_nodes[node_id];

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
                  idsComparator<Node*>);
        auto const last =
            std::unique(adjacent_nodes.begin(), adjacent_nodes.end());
        adjacent_nodes.erase(last, adjacent_nodes.end());
    }
    return nodes_connected_by_elements;
}

bool isBaseNode(Node const& node,
                std::vector<Element const*> const& elements_connected_to_node)
{
    // Check if node is connected.
    if (elements_connected_to_node.empty())
    {
        return true;
    }

    // In a mesh a node always belongs to at least one element.
    auto const e = elements_connected_to_node[0];

    auto const n_base_nodes = e->getNumberOfBaseNodes();
    auto const local_index = getNodeIDinElement(*e, &node);
    return local_index < n_base_nodes;
}

Mesh& findMeshByName(std::vector<std::unique_ptr<Mesh>> const& meshes,
                     std::string_view const name)
{
    return *BaseLib::findElementOrError(
        meshes,
        [&name](auto const& mesh)
        {
            assert(mesh != nullptr);
            return mesh->getName() == name;
        },
        [&]() { OGS_FATAL("Required mesh named {:s} not found.", name); });
}

}  // namespace MeshLib
