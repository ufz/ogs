/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>

namespace MeshLib
{

template <class ELEMENT_RULE>
TemplateElement<ELEMENT_RULE>::TemplateElement(Node* nodes[n_all_nodes], std::size_t id)
: Element(id)
{
    this->_nodes = nodes;
    this->_neighbors = new Element*[getNumberOfNeighbors()];
    std::fill(this->_neighbors, this->_neighbors + getNumberOfNeighbors(), nullptr);
    this->_content = ELEMENT_RULE::computeVolume(this->_nodes);
}

template <class ELEMENT_RULE>
TemplateElement<ELEMENT_RULE>::TemplateElement(std::array<Node*, n_all_nodes> const& nodes, std::size_t id)
: Element(id)
{
    this->_nodes = new Node*[n_all_nodes];
    std::copy(nodes.begin(), nodes.end(), this->_nodes);
    this->_neighbors = new Element*[getNumberOfNeighbors()];
    std::fill(this->_neighbors, this->_neighbors + getNumberOfNeighbors(), nullptr);
    this->_content = ELEMENT_RULE::computeVolume(this->_nodes);
}

template <class ELEMENT_RULE>
TemplateElement<ELEMENT_RULE>::TemplateElement(const TemplateElement &e)
: Element(e.getID())
{
    this->_nodes = new Node*[n_all_nodes];
    for (unsigned i=0; i<n_all_nodes; i++)
        this->_nodes[i] = e._nodes[i];
    this->_neighbors = new Element*[getNumberOfNeighbors()];
    for (unsigned i=0; i<getNumberOfNeighbors(); i++)
        this->_neighbors[i] = e._neighbors[i];
    this->_content = e.getContent();
}


namespace detail
{

template<unsigned N>
bool isEdge(unsigned const (&edge_nodes)[N], unsigned idx1, unsigned idx2)
{

    if (edge_nodes[0]==idx1 && edge_nodes[1]==idx2) return true;
    if (edge_nodes[1]==idx1 && edge_nodes[0]==idx2) return true;

    return false;
}

inline bool
isEdge(unsigned const (&/*edge_nodes*/)[1], unsigned /*idx1*/, unsigned /*idx2*/)
{
    return false;
}

} // namespace detail


template <class ELEMENT_RULE>
bool TemplateElement<ELEMENT_RULE>::isEdge(unsigned idx1, unsigned idx2) const
{
    for (unsigned i(0); i<getNumberOfEdges(); i++)
    {
        if (detail::isEdge(ELEMENT_RULE::edge_nodes[i], idx1, idx2)) return true;
    }
    return false;
}


} // MeshLib

