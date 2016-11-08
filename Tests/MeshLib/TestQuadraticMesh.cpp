/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MeshLib/MeshGenerators/QuadraticeMeshGenerator.h"
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

    auto isIn = [](MeshLib::Node const* node, std::initializer_list<std::size_t> list)
    {
        return list.end() != std::find(list.begin(), list.end(), node->getID());
    };

    for (MeshLib::Node const* node : mesh->getNodes())
    {
        if (isIn(node, {4}))
        {
            ASSERT_EQ(4u, node->getElements().size());
            ASSERT_EQ(21u, node->getConnectedNodes().size());
        }
        else if (isIn(node, {0, 2, 6, 8, 9, 10, 11, 13, 15, 18, 19, 20}))
        {
            ASSERT_EQ(1u, node->getElements().size());
            ASSERT_EQ(8u, node->getConnectedNodes().size());
        }
        else
        {
            ASSERT_EQ(2u, node->getElements().size());
            ASSERT_EQ(13u, node->getConnectedNodes().size());
        }
    }
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
        GeoLib::Polyline* ply = new GeoLib::Polyline(pnts);
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

    auto isIn = [](MeshLib::Node const* node, std::initializer_list<std::size_t> list)
    {
        return list.end() != std::find(list.begin(), list.end(), node->getID());
    };

    for (MeshLib::Node const* node : mesh->getNodes())
    {
        if (isIn(node, {4}))
        {
            ASSERT_EQ(6u, node->getElements().size());
            ASSERT_EQ(21u, node->getConnectedNodes().size());
        }
        else if (isIn(node, {0, 2, 6, 8, 9, 10, 11, 13, 15, 18, 19, 20}))
        {
            ASSERT_EQ(1u, node->getElements().size());
            ASSERT_EQ(8u, node->getConnectedNodes().size());
        }
        else if (isIn(node, {3, 5, 14, 16}))
        {
            ASSERT_EQ(3u, node->getElements().size());
            ASSERT_EQ(13u, node->getConnectedNodes().size());
        }
        else
        {
            ASSERT_EQ(2u, node->getElements().size());
            ASSERT_EQ(13u, node->getConnectedNodes().size());
        }
    }
}
