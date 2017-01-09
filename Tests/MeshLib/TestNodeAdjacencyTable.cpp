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

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"
#include "MeshLib/NodeAdjacencyTable.h"

TEST(MeshLib, CreateNodeAdjacencyTable1D)
{
    using namespace MeshLib;

    std::unique_ptr<Mesh> mesh(
        MeshGenerator::generateLineMesh(double(1), std::size_t(10)));

    NodeAdjacencyTable table(mesh->getNodes());

    // There must be as many entries as there are nodes in the mesh.
    ASSERT_EQ(mesh->getNumberOfNodes(), table.size());

    // Minimum connectivity for line mesh is 2 for boundary nodes,
    // and the maximum is 3 for interior nodes.
    std::size_t const nnodes = mesh->getNumberOfNodes();
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        std::size_t const n_connections = table.getNodeDegree(i);

        std::size_t const n_elements = mesh->getNode(i)->getNumberOfElements();
        switch (n_elements)
        {
        case 1: // a corner node has 2 adjacent nodes.
            ASSERT_EQ(2, n_connections) << " for boundary node " << i;
            break;
        case 2: // an interior node has 3 adjacent nodes.
            ASSERT_EQ(3, n_connections) << " for interior node " << i;
            break;
        default:
                FAIL() << "The regular line mesh node has unexpected number "
                    << "of elements. Node " << i << " has " << n_elements
                    << " elements.";
        }
    }
}

TEST(MeshLib, CreateNodeAdjacencyTable2D)
{
    using namespace MeshLib;

    std::unique_ptr<Mesh> mesh(MeshGenerator::generateRegularQuadMesh(
        1, 1, std::size_t(10), std::size_t(10)));

    NodeAdjacencyTable table(mesh->getNodes());

    // There must be as many entries as there are nodes in the mesh.
    ASSERT_EQ(mesh->getNumberOfNodes(), table.size());

    std::size_t const nnodes = mesh->getNumberOfNodes();
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        std::size_t const n_connections = table.getNodeDegree(i);

        std::size_t const n_elements = mesh->getNode(i)->getNumberOfElements();
        switch (n_elements)
        {
        case 1: // a corner node has 4 adjacent nodes.
                ASSERT_EQ(4, n_connections) << " for corner node " << i;
                break;
        case 2: // a boundary node has 6 adjacent nodes.
                ASSERT_EQ(6, n_connections) << " for boundary node " << i;
                break;
        case 4: // an interior node has 9 connections.
                ASSERT_EQ(9, n_connections) << " for interior node " << i;
                break;
        default:
                FAIL() << "The regular quad mesh node has unexpected number "
                    << "of elements. Node " << i << " has " << n_elements
                    << " elements.";
        }
    }
}

TEST(MeshLib, CreateNodeAdjacencyTable3D)
{
    using namespace MeshLib;

    std::unique_ptr<Mesh> mesh(MeshGenerator::generateRegularHexMesh(
                1, 1, 1, 10.0, 10.0, 10.0));
        //double(1), double(1), double(1), std::size_t(10), std::size_t(10), std::size_t(10)));

    NodeAdjacencyTable table(mesh->getNodes());

    // There must be as many entries as there are nodes in the mesh.
    ASSERT_EQ(mesh->getNumberOfNodes(), table.size());

    std::size_t const nnodes = mesh->getNumberOfNodes();
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        std::size_t const n_connections = table.getNodeDegree(i);

        std::size_t const n_elements = mesh->getNode(i)->getNumberOfElements();
        switch (n_elements)
        {
        case 1: // a corner node has 8 adjacent nodes.
                ASSERT_EQ(8, n_connections) << " for corner node " << i;
                break;
        case 2: // an edge node has 12 adjacent nodes.
                ASSERT_EQ(12, n_connections) << " for edge node " << i;
                break;
        case 3: // a boundary node has 18 adjacent nodes.
                ASSERT_EQ(18, n_connections) << " for boundary node " << i;
                break;
        case 4: // an interior node has 27 connections.
                ASSERT_EQ(27, n_connections) << " for interior node " << i;
                break;
        default:
                FAIL() << "The regular hex mesh node has unexpected number "
                    << "of elements. Node " << i << " has " << n_elements
                    << " elements.";
        }
    }
}
