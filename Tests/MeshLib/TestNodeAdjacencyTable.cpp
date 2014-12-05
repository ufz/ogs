/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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
    ASSERT_EQ(mesh->getNNodes(), table.size());

    // Minimum connectivity for line mesh is 2 for boundary nodes,
    // and the maximum is 3 for interior nodes.
    std::size_t const nnodes = mesh->getNNodes();
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        std::size_t const n_connections = table.getNodeDegree(i);
        if (i == 0 || i == mesh->getNNodes()-1)
            ASSERT_EQ(2, n_connections) << " for boundary node " << i;
        else
            ASSERT_EQ(3, n_connections) << " for interior node " << i;
    }
}

TEST(MeshLib, CreateNodeAdjacencyTable2D)
{
    using namespace MeshLib;

    std::unique_ptr<Mesh> mesh(MeshGenerator::generateRegularQuadMesh(
        double(1), double(1), std::size_t(10), std::size_t(10)));

    NodeAdjacencyTable table(mesh->getNodes());

    // There must be as many entries as there are nodes in the mesh.
    ASSERT_EQ(mesh->getNNodes(), table.size());

    std::size_t const nnodes = mesh->getNNodes();
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        std::size_t const n_connections = table.getNodeDegree(i);

        if (mesh->getNode(i)->getNElements() < 4)
        {   // boundary node has at least 4 and at most 6 adjacent nodes.
            ASSERT_LE(4, n_connections) << " for boundary node " << i;
            ASSERT_GE(6, n_connections) << " for boundary node " << i;
        }
        else
            ASSERT_EQ(9, n_connections) << " for interior node " << i;
    }
}

TEST(MeshLib, CreateNodeAdjacencyTable3D)
{
    using namespace MeshLib;

    std::unique_ptr<Mesh> mesh(MeshGenerator::generateRegularHexMesh(
                1, 1, 1, 10, 10, 10));
        //double(1), double(1), double(1), std::size_t(10), std::size_t(10), std::size_t(10)));

    NodeAdjacencyTable table(mesh->getNodes());

    // There must be as many entries as there are nodes in the mesh.
    ASSERT_EQ(mesh->getNNodes(), table.size());

    std::size_t const nnodes = mesh->getNNodes();
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        std::size_t const n_connections = table.getNodeDegree(i);

        if (mesh->getNode(i)->getNElements() < 4)
        {   // boundary node has at least 8 and at most 18 adjacent nodes.
            ASSERT_LE(8, n_connections) << " for boundary node " << i;
            ASSERT_GE(18, n_connections) << " for boundary node " << i;
        }
        else
            ASSERT_EQ(27, n_connections) << " for interior node " << i;
    }
}
