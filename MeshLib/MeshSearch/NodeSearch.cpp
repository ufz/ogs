/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NodeSearch.h"

#include <memory>
#include <set>

#include <logog/include/logog.hpp>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib {

NodeSearch::NodeSearch(const MeshLib::Mesh &mesh)
    : _mesh(mesh)
{
}

std::size_t NodeSearch::searchNodesConnectedToOnlyGivenElements(
        const std::vector<std::size_t> &elements)
{
    // Find out by how many elements a node would be removed.
    //
    // Note: If there are only few elements to be removed, using a different
    // algorithm might be more memory efficient.
    std::vector<std::size_t> node_marked_counts(_mesh.getNumberOfNodes(), 0);

    for(std::size_t eid : elements)
    {
        auto* e = _mesh.getElement(eid);
        for (unsigned i=0; i<e->getNumberOfNodes(); i++) {
            node_marked_counts[e->getNodeIndex(i)]++;
        }
    }


    // Push back nodes which counts are equal to number of connected elements to
    // that node.
    std::vector<std::size_t> connected_nodes;
    for (std::size_t i=0; i<node_marked_counts.size(); i++)
    {
        if (node_marked_counts[i] == _mesh.getNode(i)->getElements().size())
            connected_nodes.push_back(i);
    }

    this->updateUnion(connected_nodes);
    return connected_nodes.size();
}

std::size_t NodeSearch::searchUnused()
{
    const std::size_t nNodes (_mesh.getNumberOfNodes());
    const std::vector<MeshLib::Node*> &nodes (_mesh.getNodes());
    std::vector<std::size_t> del_node_idx;

    for (unsigned i=0; i<nNodes; ++i)
        if (nodes[i]->getNumberOfElements() == 0)
            del_node_idx.push_back(i);

    this->updateUnion(del_node_idx);
    return del_node_idx.size();
}

std::size_t NodeSearch::searchBoundaryNodes()
{
    std::vector<std::size_t> vec_boundary_nodes;
    if (_mesh.getDimension() == 1)
    {
        for (MeshLib::Node const* n : _mesh.getNodes())
            if (n->getElements().size() == 1)
                vec_boundary_nodes.push_back(n->getID());
    }
    else if (_mesh.getDimension() == 2)
    {
        for (MeshLib::Element const* elem : _mesh.getElements())
        {
            if (elem->getDimension() < _mesh.getDimension())
                continue;
            if (!elem->isBoundaryElement())
                continue;

            std::size_t const n_edges (elem->getNumberOfEdges());
            for (std::size_t i=0; i<n_edges; ++i)
            {
                if (elem->getNeighbor(i) != nullptr)
                    continue;
                std::unique_ptr<MeshLib::Element const> edge(elem->getEdge(i));
                for (unsigned j=0; j<edge->getNumberOfNodes(); j++)
                    vec_boundary_nodes.push_back(edge->getNode(j)->getID());
            }
        }
    }
    else
    {
        for (MeshLib::Element const* elem : _mesh.getElements())
        {
            if (elem->getDimension() < _mesh.getDimension())
                continue;
            if (!elem->isBoundaryElement())
                continue;

            std::size_t const n_faces (elem->getNumberOfFaces());
            for (std::size_t i=0; i<n_faces; ++i)
            {
                if (elem->getNeighbor(i) != nullptr)
                    continue;
                std::unique_ptr<MeshLib::Element const> face(elem->getFace(i));
                for (unsigned j=0; j<face->getNumberOfNodes(); j++)
                    vec_boundary_nodes.push_back(face->getNode(j)->getID());
            }
        }
    }
    std::sort(vec_boundary_nodes.begin(), vec_boundary_nodes.end());
    vec_boundary_nodes.erase(std::unique(vec_boundary_nodes.begin(), vec_boundary_nodes.end()), vec_boundary_nodes.end());


    this->updateUnion(vec_boundary_nodes);
    return vec_boundary_nodes.size();
}

void NodeSearch::updateUnion(const std::vector<std::size_t> &vec)
{
    std::vector<std::size_t> vec_temp(vec.size() + _marked_nodes.size());
    auto it = std::set_union(vec.begin(), vec.end(), _marked_nodes.begin(), _marked_nodes.end(), vec_temp.begin());
    vec_temp.resize(it - vec_temp.begin());
    _marked_nodes.assign(vec_temp.begin(), vec_temp.end());
}

std::vector<Node*> getUniqueNodes(std::vector<Element*> const& elements)
{
    std::set<Node*> nodes_set;
    for (auto e : elements)
    {
        Node* const* nodes = e->getNodes();
        unsigned const nnodes = e->getNumberOfNodes();
        nodes_set.insert(nodes, nodes + nnodes);
    }

    std::vector<Node*> nodes;
    nodes.reserve(nodes_set.size());

    std::move(nodes_set.cbegin(), nodes_set.cend(),
        std::back_inserter(nodes));

    return nodes;
}

} // end namespace MeshLib
