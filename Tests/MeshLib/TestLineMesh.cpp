/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <ctime>

#include "MeshLib/Elements/Line.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

class MeshLibLineMesh : public ::testing::Test
{
public:
    MeshLibLineMesh()
    {
        mesh = MeshToolsLib::MeshGenerator::generateLineMesh(extent, mesh_size);
    }

    ~MeshLibLineMesh() override { delete mesh; }
    static std::size_t const mesh_size = 9;
    double extent = 1.0;
    MeshLib::Mesh const* mesh{nullptr};
};
std::size_t const MeshLibLineMesh::mesh_size;

TEST_F(MeshLibLineMesh, Construction)
{
    ASSERT_TRUE(mesh != nullptr);

    // There are mesh_size elements in the mesh.
    ASSERT_EQ(mesh_size, mesh->getNumberOfElements());

    // There are mesh_size+1 nodes in the mesh.
    ASSERT_EQ(mesh_size + 1, mesh->getNumberOfNodes());

    // All elements have maximum two neighbors.
    std::vector<MeshLib::Element*> const& elements = mesh->getElements();
    for (auto e : elements)
    {
        ASSERT_EQ(2u, e->getNumberOfNeighbors());
    }

    auto const [min, max] = minMaxEdgeLength(mesh->getElements());
    ASSERT_NEAR(extent / mesh_size, min,
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(extent / mesh_size, max,
                std::numeric_limits<double>::epsilon());
}

TEST_F(MeshLibLineMesh, ElementNeigbors)
{
    auto count_neighbors = [](MeshLib::Element const* const e)
    {
        unsigned count = 0;
        for (int i = 0; i < 2; i++)
        {
            if (e->getNeighbor(i) != nullptr)
            {
                count++;
            }
        }
        return count;
    };

    std::vector<MeshLib::Element*> const& elements = mesh->getElements();

    // Each element has two neighbors n-1 and n+1 except the first and the last
    // elements.
    ASSERT_EQ(1u, count_neighbors(elements.front()));
    ASSERT_EQ(1u, count_neighbors(elements.back()));
    ASSERT_TRUE(areNeighbors(elements.front(), elements[1]));
    ASSERT_TRUE(areNeighbors(elements.back(), elements[elements.size() - 2]));

    for (std::size_t i = 1; i < elements.size() - 1; ++i)
    {
        ASSERT_EQ(2u, count_neighbors(elements[i]));
        ASSERT_TRUE(areNeighbors(elements[i], elements[i - 1]));
        ASSERT_TRUE(areNeighbors(elements[i], elements[i + 1]));
    }
}

TEST_F(MeshLibLineMesh, ElementToNodeConnectivity)
{
    std::vector<MeshLib::Element*> const& elements = mesh->getElements();
    // Each element n of the line mesh is formed by nodes n and n+1.
    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        // An element consists of two nodes n and n+1
        ASSERT_EQ(2u, elements[i]->getNumberOfBaseNodes());
        ASSERT_EQ(i, elements[i]->getNode(0)->getID());
        ASSERT_EQ(i + 1, elements[i]->getNode(1)->getID());
    }
}

TEST_F(MeshLibLineMesh, NodeToElementConnectivity)
{
    std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
    std::vector<MeshLib::Element*> const& elements = mesh->getElements();

    // Each node n of the line mesh is connected to two elements n-1 and n,
    // except the first and last nodes.
    ASSERT_EQ(1u, mesh->getElementsConnectedToNode(*nodes.front()).size());
    ASSERT_EQ(elements[0], mesh->getElementsConnectedToNode(*nodes.front())[0]);

    ASSERT_EQ(1u, mesh->getElementsConnectedToNode(*nodes.back()).size());
    ASSERT_EQ(elements.back(),
              mesh->getElementsConnectedToNode(*nodes.back())[0]);

    for (std::size_t i = 1; i < nodes.size() - 1; ++i)
    {
        ASSERT_EQ(2u, mesh->getElementsConnectedToNode(*nodes[i]).size());
        ASSERT_EQ(elements[i - 1],
                  mesh->getElementsConnectedToNode(*nodes[i])[0]);
        ASSERT_EQ(elements[i], mesh->getElementsConnectedToNode(*nodes[i])[1]);
    }
}
