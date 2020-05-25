/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    this->nodes_ = nodes;
    this->neighbors_ = new Element*[getNumberOfNeighbors()];
    std::fill(this->neighbors_, this->neighbors_ + getNumberOfNeighbors(), nullptr);
    this->content_ = ELEMENT_RULE::computeVolume(this->nodes_);
}

template <class ELEMENT_RULE>
TemplateElement<ELEMENT_RULE>::TemplateElement(std::array<Node*, n_all_nodes> const& nodes, std::size_t id)
: Element(id)
{
    this->nodes_ = new Node*[n_all_nodes];
    std::copy(nodes.begin(), nodes.end(), this->nodes_);
    this->neighbors_ = new Element*[getNumberOfNeighbors()];
    std::fill(this->neighbors_, this->neighbors_ + getNumberOfNeighbors(), nullptr);
    this->content_ = ELEMENT_RULE::computeVolume(this->nodes_);
}

template <class ELEMENT_RULE>
TemplateElement<ELEMENT_RULE>::TemplateElement(const TemplateElement &e)
: Element(e.getID())
{
    this->nodes_ = new Node*[n_all_nodes];
    for (unsigned i = 0; i < n_all_nodes; i++)
    {
        this->nodes_[i] = e.nodes_[i];
    }
    this->neighbors_ = new Element*[getNumberOfNeighbors()];
    for (unsigned i = 0; i < getNumberOfNeighbors(); i++)
    {
        this->neighbors_[i] = e.neighbors_[i];
    }
    this->content_ = e.getContent();
}


namespace details
{

template<unsigned N>
bool isEdge(unsigned const (&edge_nodes)[N], unsigned idx1, unsigned idx2)
{
    if (edge_nodes[0] == idx1 && edge_nodes[1] == idx2)
    {
        return true;
    }
    if (edge_nodes[1] == idx1 && edge_nodes[0] == idx2)
    {
        return true;
    }

    return false;
}

inline bool
isEdge(unsigned const (&/*edge_nodes*/)[1], unsigned /*idx1*/, unsigned /*idx2*/)
{
    return false;
}

} // namespace details


template <class ELEMENT_RULE>
bool TemplateElement<ELEMENT_RULE>::isEdge(unsigned idx1, unsigned idx2) const
{
    for (unsigned i(0); i<getNumberOfEdges(); i++)
    {
        if (details::isEdge(ELEMENT_RULE::edge_nodes[i], idx1, idx2))
        {
            return true;
        }
    }
    return false;
}

}  // namespace MeshLib
