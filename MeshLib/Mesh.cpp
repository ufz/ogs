/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief Implementation of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Mesh.h"

#include <memory>

#include "BaseLib/RunTime.h"

#include "Elements/Element.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"

namespace MeshLib
{

Mesh::Mesh(const std::string &name,
           const std::vector<Node*> &nodes,
           const std::vector<Element*> &elements,
           Properties const& properties,
           const std::size_t n_base_nodes)
    : _id(_counter_value-1), _mesh_dimension(0),
      _edge_length(std::numeric_limits<double>::max(), 0),
      _node_distance(std::numeric_limits<double>::max(), 0),
      _name(name), _nodes(nodes), _elements(elements),
      _n_base_nodes(n_base_nodes),
      _properties(properties)
{
    assert(n_base_nodes <= nodes.size());
    this->resetNodeIDs();
    this->resetElementIDs();
    if ((n_base_nodes==0 && hasNonlinearElement()) || isNonlinear())
        this->checkNonlinearNodeIDs();
    this->setDimension();
    this->setElementsConnectedToNodes();
    //this->setNodesConnectedByEdges();
    this->setNodesConnectedByElements();
    this->setElementNeighbors();

    this->calcEdgeLengthRange();
}

Mesh::Mesh(const Mesh &mesh)
    : _id(_counter_value-1), _mesh_dimension(mesh.getDimension()),
      _edge_length(mesh._edge_length.first, mesh._edge_length.second),
      _node_distance(mesh._node_distance.first, mesh._node_distance.second),
      _name(mesh.getName()), _nodes(mesh.getNumberOfNodes()), _elements(mesh.getNumberOfElements()),
      _n_base_nodes(mesh.getNumberOfBaseNodes()),
      _properties(mesh._properties)
{
    const std::vector<Node*>& nodes (mesh.getNodes());
    const std::size_t nNodes (nodes.size());
    for (unsigned i=0; i<nNodes; ++i)
        _nodes[i] = new Node(*nodes[i]);

    const std::vector<Element*>& elements (mesh.getElements());
    const std::size_t nElements (elements.size());
    for (unsigned i=0; i<nElements; ++i)
    {
        const std::size_t nElemNodes = elements[i]->getNumberOfNodes();
        _elements[i] = elements[i]->clone();
        for (unsigned j=0; j<nElemNodes; ++j)
            _elements[i]->_nodes[j] = _nodes[elements[i]->getNode(j)->getID()];
    }

    if (_mesh_dimension==0) this->setDimension();
    this->setElementsConnectedToNodes();
    //this->setNodesConnectedByEdges();
    //this->setNodesConnectedByElements();
    this->setElementNeighbors();
}

Mesh::~Mesh()
{
    const std::size_t nElements (_elements.size());
    for (std::size_t i=0; i<nElements; ++i)
        delete _elements[i];

    const std::size_t nNodes (_nodes.size());
    for (std::size_t i=0; i<nNodes; ++i)
        delete _nodes[i];
}

void Mesh::addNode(Node* node)
{
    _nodes.push_back(node);
}

void Mesh::addElement(Element* elem)
{
    _elements.push_back(elem);

    // add element information to nodes
    unsigned nNodes (elem->getNumberOfNodes());
    for (unsigned i=0; i<nNodes; ++i)
        elem->_nodes[i]->addElement(elem);
}

void Mesh::resetNodeIDs()
{
    const std::size_t nNodes (this->_nodes.size());
    for (unsigned i=0; i<nNodes; ++i)
        _nodes[i]->setID(i);

    if (_n_base_nodes==0)
    {
        unsigned max_basenode_ID = 0;
        for (Element const* e : _elements)
            for (unsigned i=0; i<e->getNumberOfBaseNodes(); i++)
                max_basenode_ID = std::max(max_basenode_ID, e->getNodeIndex(i));
        _n_base_nodes = max_basenode_ID + 1;
    }
}

void Mesh::resetElementIDs()
{
    const std::size_t nElements (this->_elements.size());
    for (unsigned i=0; i<nElements; ++i)
        _elements[i]->setID(i);
}

void Mesh::setDimension()
{
    const std::size_t nElements (_elements.size());
    for (unsigned i=0; i<nElements; ++i)
        if (_elements[i]->getDimension() > _mesh_dimension)
            _mesh_dimension = _elements[i]->getDimension();
}

void Mesh::setElementsConnectedToNodes()
{
    for (auto e = _elements.begin(); e != _elements.end(); ++e)
    {
        const unsigned nNodes ((*e)->getNumberOfNodes());
        for (unsigned j=0; j<nNodes; ++j)
            (*e)->_nodes[j]->addElement(*e);
    }
}

void Mesh::resetElementsConnectedToNodes()
{
    for (auto node = _nodes.begin(); node != _nodes.end(); ++node)
        if (*node)
            (*node)->clearElements();
    this->setElementsConnectedToNodes();
}

void Mesh::calcEdgeLengthRange()
{
    this->_edge_length.first  = std::numeric_limits<double>::max();
    this->_edge_length.second = 0;
    double min_length(0);
    double max_length(0);
    const std::size_t nElems (this->getNumberOfElements());
    for (std::size_t i=0; i<nElems; ++i)
    {
        _elements[i]->computeSqrEdgeLengthRange(min_length, max_length);
        this->_edge_length.first  = std::min(this->_edge_length.first,  min_length);
        this->_edge_length.second = std::max(this->_edge_length.second, max_length);
    }
    this->_edge_length.first  = sqrt(this->_edge_length.first);
    this->_edge_length.second = sqrt(this->_edge_length.second);
}

void Mesh::setElementNeighbors()
{
    std::vector<Element*> neighbors;
    for (auto it = _elements.begin(); it != _elements.end(); ++it)
    {
        // create vector with all elements connected to current element (includes lots of doubles!)
        Element *const element = *it;

        const std::size_t nNodes (element->getNumberOfBaseNodes());
        for (unsigned n(0); n<nNodes; ++n)
        {
            std::vector<Element*> const& conn_elems ((element->getNode(n)->getElements()));
            neighbors.insert(neighbors.end(), conn_elems.begin(), conn_elems.end());
        }
        std::sort(neighbors.begin(), neighbors.end());
        auto const neighbors_new_end = std::unique(neighbors.begin(), neighbors.end());

        for (auto neighbor = neighbors.begin(); neighbor != neighbors_new_end; ++neighbor)
        {
            boost::optional<unsigned> const opposite_face_id = element->addNeighbor(*neighbor);
            if (opposite_face_id)
            {
                (*neighbor)->setNeighbor(element, *opposite_face_id);
            }
        }
        neighbors.clear();
    }
}

void Mesh::setNodesConnectedByEdges()
{
    const std::size_t nNodes (this->_nodes.size());
    for (unsigned i=0; i<nNodes; ++i)
    {
        MeshLib::Node* node (_nodes[i]);
        std::vector<MeshLib::Node*> conn_set;
        const std::vector<MeshLib::Element*> &conn_elems (node->getElements());
        const std::size_t nConnElems (conn_elems.size());
        for (unsigned j=0; j<nConnElems; ++j)
        {
            MeshLib::Element* conn_ele = conn_elems[j];
            const unsigned idx (conn_ele->getNodeIDinElement(node));
            const unsigned nElemNodes (conn_ele->getNumberOfBaseNodes());
            for (unsigned k(0); k<nElemNodes; ++k)
            {
                MeshLib::Node const* node_k = conn_ele->getNode(k);
                bool is_in_vector (false);
                const std::size_t nConnNodes (conn_set.size());
                for (unsigned l(0); l<nConnNodes; ++l)
                    if (node_k == conn_set[l])
                        is_in_vector = true;
                if (is_in_vector) continue;

                if (conn_ele->getNumberOfBaseNodes() == conn_ele->getNumberOfNodes())
                {
                    if (conn_ele->isEdge(idx, k))
                        conn_set.push_back(const_cast<MeshLib::Node*>(node_k));
                }
                else
                {
                    for (unsigned l=0; l<conn_ele->getNumberOfEdges(); l++)
                    {
                        std::unique_ptr<Element const> edge(conn_ele->getEdge(l));
                        unsigned match = 0;
                        for (unsigned m=0; m<edge->getNumberOfBaseNodes(); m++)
                        {
                            auto edge_node = edge->getNode(m);
                            if (edge_node == node || edge_node == node_k)
                                match++;
                        }
                        if (match != 2)
                            continue;
                        conn_set.push_back(const_cast<MeshLib::Node*>(node_k));
                        for (unsigned m=edge->getNumberOfBaseNodes(); m<edge->getNumberOfNodes(); m++)
                            conn_set.push_back(const_cast<MeshLib::Node*>(edge->getNode(m)));
                        break;
                    }
                }

            }
        }
        node->setConnectedNodes(conn_set);
    }
}

void Mesh::setNodesConnectedByElements()
{
    // Allocate temporary space for adjacent nodes.
    std::vector<Node*> adjacent_nodes;
    for (Node* const node : _nodes)
    {
        adjacent_nodes.clear();

        // Get all elements, to which this node is connected.
        std::vector<Element*> const& conn_elems = node->getElements();

        // And collect all elements' nodes.
        for (Element const* const element : conn_elems)
        {
            Node* const* const single_elem_nodes = element->getNodes();
            std::size_t const nnodes = element->getNumberOfNodes();
            for (std::size_t n = 0; n < nnodes; n++)
                adjacent_nodes.push_back(single_elem_nodes[n]);
        }

        // Make nodes unique and sorted by their ids.
        // This relies on the node's id being equivalent to it's address.
        std::sort(adjacent_nodes.begin(), adjacent_nodes.end(),
            [](Node* a, Node* b) { return a->getID() < b->getID(); });
        auto const last = std::unique(adjacent_nodes.begin(), adjacent_nodes.end());
        adjacent_nodes.erase(last, adjacent_nodes.end());

        node->setConnectedNodes(adjacent_nodes);
    }
}

void Mesh::checkNonlinearNodeIDs() const
{
    for (MeshLib::Element const* e : _elements)
    {
        for (unsigned i=e->getNumberOfBaseNodes(); i<e->getNumberOfNodes(); i++)
        {
            if (e->getNodeIndex(i) >= getNumberOfBaseNodes())
                continue;

            ERR("Found a nonlinear node whose ID (%d) is smaller than the "
                "number of base node IDs (%d)."
                "Some functions may not work properly.",
                e->getNodeIndex(i), getNumberOfBaseNodes());
            return;
        }
    }
}

bool Mesh::hasNonlinearElement() const
{
    return std::any_of(std::begin(_elements), std::end(_elements),
        [](Element const* const e) {
            return e->getNumberOfNodes() != e->getNumberOfBaseNodes();
        });
}

void scaleMeshPropertyVector(MeshLib::Mesh & mesh,
                             std::string const& property_name,
                             double factor)
{
    for (auto& v :
         *mesh.getProperties().getPropertyVector<double>(property_name))
        v *= factor;
}

}
