/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include "QuadraticeMeshGenerator.h"

#include "BaseLib/makeVectorUnique.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"

namespace MeshLib
{
namespace
{

struct Edge
{
    Edge(std::size_t id1, std::size_t id2)
        : _id1(std::min(id1, id2)), _id2(std::max(id1, id2))
    {}

    bool operator==(Edge const& r) const
    {
        return (_id1 == r._id1 && _id2 == r._id2);
    }

    std::size_t _id1;
    std::size_t _id2;

    std::size_t _edge_id = 0;
};

bool operator< (Edge const& l, Edge const& r)
{
    return (l._id1 != r._id1) ? (l._id1 < r._id1) : l._id2 < r._id2;
}

template <typename T_ELEMENT>
T_ELEMENT* createQuadraticElement(MeshLib::Element const* e, std::vector<Edge> const& vec_edges,
                  std::vector<MeshLib::Node*> const& vec_new_nodes, const std::size_t n_mesh_base_nodes)
{
    auto const n_all_nodes = T_ELEMENT::n_all_nodes;
    auto const n_base_nodes = T_ELEMENT::n_base_nodes;
    MeshLib::Node** nodes = new MeshLib::Node*[n_all_nodes];
    for (unsigned i=0; i<e->getNumberOfBaseNodes(); i++)
        nodes[i] = const_cast<MeshLib::Node*>(vec_new_nodes[e->getNode(i)->getID()]);
    for (unsigned i=0; i<e->getNumberOfEdges(); i++)
    {
        auto itr = std::find(vec_edges.begin(), vec_edges.end(),
                             Edge(e->getEdgeNode(i, 0)->getID(),
                                  e->getEdgeNode(i, 1)->getID()));
        assert(itr != vec_edges.end());
        nodes[n_base_nodes + i] =
            vec_new_nodes[n_mesh_base_nodes + itr->_edge_id];
    }
    return new T_ELEMENT(nodes);
}

} // no named namespace

std::unique_ptr<Mesh> createQuadraticOrderMesh(Mesh const& org_mesh)
{
    std::vector<MeshLib::Node*> vec_new_nodes = MeshLib::copyNodeVector(org_mesh.getNodes());

    // collect edges
    std::vector<Edge> vec_edges;
    for (MeshLib::Element const* e : org_mesh.getElements())
    {
        for (unsigned i=0; i<e->getNumberOfEdges(); i++)
        {
            auto node0 = e->getEdgeNode(i, 0);
            auto node1 = e->getEdgeNode(i, 1);
            vec_edges.emplace_back(node0->getID(), node1->getID());
        }
    }
    BaseLib::makeVectorUnique(vec_edges);
    for (std::size_t i=0; i<vec_edges.size(); i++)
         vec_edges[i]._edge_id = i;
    INFO("Found %d edges in the mesh", vec_edges.size());

    // create mid-point nodes
    double coords[3];
    for (Edge const& edge : vec_edges)
    {
        auto const& node0 = *vec_new_nodes[edge._id1];
        auto const& node1 = *vec_new_nodes[edge._id2];
        for (unsigned i=0; i<3; i++)
            coords[i] = (node0[i] + node1[i]) * 0.5;
        vec_new_nodes.push_back(new MeshLib::Node(coords));
    }

    // create new elements with the quadratic nodes
    std::vector<MeshLib::Element*> vec_new_eles;
    for (MeshLib::Element const* e : org_mesh.getElements())
    {
        if (e->getCellType() == MeshLib::CellType::LINE2)
        {
            vec_new_eles.push_back(createQuadraticElement<MeshLib::Line3>(
                e, vec_edges, vec_new_nodes, org_mesh.getNumberOfNodes()));
        }
        else if (e->getCellType() == MeshLib::CellType::QUAD4)
        {
            vec_new_eles.push_back(createQuadraticElement<MeshLib::Quad8>(
                e, vec_edges, vec_new_nodes, org_mesh.getNumberOfNodes()));
        }
        else
        {
            OGS_FATAL("Mesh element type %s is not supported", MeshLib::CellType2String(e->getCellType()).c_str());
        }
    }

    std::unique_ptr<MeshLib::Mesh> new_mesh(
        new MeshLib::Mesh(org_mesh.getName(), vec_new_nodes, vec_new_eles,
                          org_mesh.getProperties().excludeCopyProperties(
                              std::vector<MeshLib::MeshItemType>(
                                  1, MeshLib::MeshItemType::Node)),
                          org_mesh.getNumberOfNodes()));
    return new_mesh;
}

} // namespace MeshLib

