/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>

namespace MeshLib
{
template <class ELEMENT_RULE>
TemplateElement<ELEMENT_RULE>::TemplateElement(Node* nodes[n_all_nodes],
                                               std::size_t id)
    : Element(id)
{
    std::copy_n(nodes, n_all_nodes, std::begin(_nodes));
    delete[] nodes;
    this->_neighbors = new Element*[ELEMENT_RULE::n_neighbors];
    std::fill(this->_neighbors, this->_neighbors + ELEMENT_RULE::n_neighbors,
              nullptr);

    this->space_dimension_ = ELEMENT_RULE::dimension;
}

template <class ELEMENT_RULE>
TemplateElement<ELEMENT_RULE>::TemplateElement(
    std::array<Node*, n_all_nodes> const& nodes, std::size_t id)
    : Element(id), _nodes{nodes}
{
    this->_neighbors = new Element*[ELEMENT_RULE::n_neighbors];
    std::fill(this->_neighbors, this->_neighbors + ELEMENT_RULE::n_neighbors,
              nullptr);

    this->space_dimension_ = ELEMENT_RULE::dimension;
}

template <class ELEMENT_RULE>
TemplateElement<ELEMENT_RULE>::TemplateElement(
    TemplateElement<ELEMENT_RULE> const& e)
    : Element(e.getID()), _nodes{e._nodes}
{
    this->_neighbors = new Element*[ELEMENT_RULE::n_neighbors];
    for (unsigned i = 0; i < ELEMENT_RULE::n_neighbors; i++)
    {
        this->_neighbors[i] = e._neighbors[i];
    }

    this->space_dimension_ = e.space_dimension_;
}

template <class ELEMENT_RULE>
double TemplateElement<ELEMENT_RULE>::getContent() const
{
    return ELEMENT_RULE::computeVolume(_nodes.data());
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

template <class ELEMENT_RULE>
const Node* TemplateElement<ELEMENT_RULE>::getNode(unsigned const idx) const
{
#ifndef NDEBUG
    if (idx >= getNumberOfNodes())
    {
        ERR("Error in MeshLib::TemplateElement::getNode() - Index {:d} in {:s}",
            idx, MeshElemType2String(getGeomType()));
        return nullptr;
    }
#endif
    return _nodes[idx];
}

template <class ELEMENT_RULE>
Node* TemplateElement<ELEMENT_RULE>::getNode(unsigned const idx)
{
#ifndef NDEBUG
    if (idx >= getNumberOfNodes())
    {
        ERR("Error in MeshLib::TemplateElement::getNode() - Index {:d} in {:s}",
            idx, MeshElemType2String(getGeomType()));
        return nullptr;
    }
#endif
    return _nodes[idx];
}

template <class ELEMENT_RULE>
void TemplateElement<ELEMENT_RULE>::setNode(unsigned idx, Node* node)
{
#ifndef NDEBUG
    if (idx < getNumberOfNodes())
#endif
    {
        _nodes[idx] = node;
    }
}

}  // namespace MeshLib
