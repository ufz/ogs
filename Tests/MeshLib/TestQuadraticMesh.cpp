/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "gtest/gtest.h"

#include <memory>

#include "GeoLib/Polyline.h"
#include "GeoLib/PolylineVec.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshGenerators/QuadraticMeshGenerator.h"
#include "MeshLib/Node.h"
#include "MeshGeoToolsLib/AppendLinesAlongPolyline.h"

TEST(MeshLib, QuadraticOrderMesh_Line)
{
    using namespace MeshLib;

    std::unique_ptr<Mesh> linear_mesh(MeshGenerator::generateLineMesh(
        1, std::size_t(2)));
    std::unique_ptr<Mesh> mesh(createQuadraticOrderMesh(*linear_mesh.get()));
    ASSERT_EQ(5u, mesh->getNumberOfNodes());
    ASSERT_EQ(3u, mesh->getNumberOfBaseNodes());
    ASSERT_EQ(2u, mesh->getNumberOfElements());

    for (MeshLib::Element const* e : mesh->getElements())
    {
        ASSERT_EQ(MeshLib::CellType::LINE3, e->getCellType());
        ASSERT_EQ(2u, e->getNumberOfBaseNodes());
        ASSERT_EQ(3u, e->getNumberOfNodes());

        for (unsigned i=0; i<e->getNumberOfBaseNodes(); i++)
            ASSERT_TRUE(mesh->isBaseNode(e->getNodeIndex(i)));
        for (unsigned i=e->getNumberOfBaseNodes(); i<e->getNumberOfNodes(); i++)
            ASSERT_FALSE(mesh->isBaseNode(e->getNodeIndex(i)));
    }

    for (MeshLib::Node const* node : mesh->getNodes())
    {
        if (node->getID() == 1)
        {
            ASSERT_EQ(2u, node->getElements().size());
            ASSERT_EQ(5u, node->getConnectedNodes().size());
        }
        else
        {
            ASSERT_EQ(1u, node->getElements().size());
            ASSERT_EQ(3u, node->getConnectedNodes().size());
        }
    }
}

TEST(MeshLib, QuadraticOrderMesh_Quad)
{
    using namespace MeshLib;

    std::unique_ptr<Mesh> linear_mesh(MeshGenerator::generateRegularQuadMesh(
        1, 1, std::size_t(2), std::size_t(2)));
    std::unique_ptr<Mesh> mesh(createQuadraticOrderMesh(*linear_mesh.get()));
    ASSERT_EQ(21u, mesh->getNumberOfNodes());
    ASSERT_EQ(9u, mesh->getNumberOfBaseNodes());
    ASSERT_EQ(4u, mesh->getNumberOfElements());

    for (MeshLib::Element const* e : mesh->getElements())
    {
        ASSERT_EQ(MeshLib::CellType::QUAD8, e->getCellType());
        ASSERT_EQ(4u, e->getNumberOfBaseNodes());
        ASSERT_EQ(8u, e->getNumberOfNodes());

        for (unsigned i=0; i<e->getNumberOfBaseNodes(); i++)
            ASSERT_TRUE(mesh->isBaseNode(e->getNodeIndex(i)));
        for (unsigned i=e->getNumberOfBaseNodes(); i<e->getNumberOfNodes(); i++)
            ASSERT_FALSE(mesh->isBaseNode(e->getNodeIndex(i)));
    }

    auto const& mesh_nodes = mesh->getNodes();

    // Count nodes shared by four elements and also connected to all 21 other
    // nodes.
    ASSERT_EQ(1, std::count_if(mesh_nodes.begin(), mesh_nodes.end(),
                               [](Node* const n) {
                                   return (n->getElements().size() == 4) &&
                                          (n->getConnectedNodes().size() == 21);
                               }));

    // Count nodes belonging to one element and also connected to all 8 other
    // nodes of that corner element.
    ASSERT_EQ(12, std::count_if(mesh_nodes.begin(), mesh_nodes.end(),
                                [](Node* const n) {
                                    return (n->getElements().size() == 1) &&
                                           (n->getConnectedNodes().size() == 8);
                                }));

    // Count nodes shared by two elements and also connected to the 13 other
    // nodes of the two elements.
    ASSERT_EQ(8, std::count_if(mesh_nodes.begin(), mesh_nodes.end(),
                               [](Node* const n) {
                                   return (n->getElements().size() == 2) &&
                                          (n->getConnectedNodes().size() == 13);
                               }));
}

TEST(MeshLib, QuadraticOrderMesh_LineQuad)
{
    using namespace MeshLib;

    std::unique_ptr<Mesh> linear_quad_mesh(MeshGenerator::generateRegularQuadMesh(
        1, 1, std::size_t(2), std::size_t(2)));

    std::unique_ptr<Mesh> linear_mesh;
    {
        std::vector<GeoLib::Point*> pnts;
        pnts.push_back(new GeoLib::Point(0,0.5,0,0));
        pnts.push_back(new GeoLib::Point(1,0.5,0,1));
        auto* ply = new GeoLib::Polyline(pnts);
        ply->addPoint(0);
        ply->addPoint(1);
        std::unique_ptr<std::vector<GeoLib::Polyline*>> plys(new std::vector<GeoLib::Polyline*>());
        plys->push_back(ply);
        GeoLib::PolylineVec ply_vec("", std::move(plys));

        linear_mesh = MeshGeoToolsLib::appendLinesAlongPolylines(*linear_quad_mesh.get(), ply_vec);

        for (auto p : pnts)
            delete p;
    }
    ASSERT_EQ(6u, linear_mesh->getNumberOfElements());

    std::unique_ptr<Mesh> mesh(createQuadraticOrderMesh(*linear_mesh.get()));
    ASSERT_EQ(21u, mesh->getNumberOfNodes());
    ASSERT_EQ(9u, mesh->getNumberOfBaseNodes());
    ASSERT_EQ(6u, mesh->getNumberOfElements());

    for (MeshLib::Element const* e : mesh->getElements())
    {
        if (e->getID() < 4)
        {
            ASSERT_EQ(MeshLib::CellType::QUAD8, e->getCellType());
            ASSERT_EQ(4u, e->getNumberOfBaseNodes());
            ASSERT_EQ(8u, e->getNumberOfNodes());
        }
        else
        {
            ASSERT_EQ(MeshLib::CellType::LINE3, e->getCellType());
            ASSERT_EQ(2u, e->getNumberOfBaseNodes());
            ASSERT_EQ(3u, e->getNumberOfNodes());
        }

        for (unsigned i=0; i<e->getNumberOfBaseNodes(); i++)
            ASSERT_TRUE(mesh->isBaseNode(e->getNodeIndex(i)));
        for (unsigned i=e->getNumberOfBaseNodes(); i<e->getNumberOfNodes(); i++)
            ASSERT_FALSE(mesh->isBaseNode(e->getNodeIndex(i)));
    }

    auto const& mesh_nodes = mesh->getNodes();

    // Count nodes shared by six elements and also connected to all 21 other
    // nodes.
    ASSERT_EQ(1, std::count_if(mesh_nodes.begin(), mesh_nodes.end(),
                               [](Node* const n) {
                                   return (n->getElements().size() == 6) &&
                                          (n->getConnectedNodes().size() == 21);
                               }));

    // Count nodes belonging to one element and also connected to all 8 other
    // nodes of that corner element.
    ASSERT_EQ(12, std::count_if(mesh_nodes.begin(), mesh_nodes.end(),
                                [](Node* const n) {
                                    return (n->getElements().size() == 1) &&
                                           (n->getConnectedNodes().size() == 8);
                                }));

    // Count nodes shared by three elements (quads and the line) and also
    // connected to the 13 other nodes of the two elements.
    ASSERT_EQ(4, std::count_if(mesh_nodes.begin(), mesh_nodes.end(),
                               [](Node* const n) {
                                   return (n->getElements().size() == 3) &&
                                          (n->getConnectedNodes().size() == 13);
                               }));

    // Count nodes shared by two elements (quads) and also connected to the 13
    // other nodes of the two elements.
    ASSERT_EQ(4, std::count_if(mesh_nodes.begin(), mesh_nodes.end(),
                               [](Node* const n) {
                                   return (n->getElements().size() == 2) &&
                                          (n->getConnectedNodes().size() == 13);
                               }));
}
