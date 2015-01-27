/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLIB_COMPONENTSELECTOR_H_
#define MESHLIB_COMPONENTSELECTOR_H_

namespace MeshLib
{


/// Create a vector of unique nodes used by given elements.
inline
std::vector<MeshLib::Node*>
selectNodes(std::vector<MeshLib::Element*> const& elements)
{
    std::set<MeshLib::Node*> nodes_set;
    for (auto e : elements)
    {
        MeshLib::Node* const* nodes = e->getNodes();
        unsigned const nnodes = e->getNNodes();
        nodes_set.insert(nodes, nodes + nnodes);
    }

    std::vector<MeshLib::Node*> nodes;
    nodes.reserve(nodes_set.size());

    std::move(nodes_set.cbegin(), nodes_set.cend(),
        std::back_inserter(nodes));

    return nodes;
}

}   // namespace MeshLib

#endif  // MESHLIB_COMPONENTSELECTOR_H_
