/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Line.h"

/// Creates a mesh consisiting of two triangles sharing a line element on their
/// common edge.
class MeshLibTriLineMesh : public ::testing::Test
{
    public:
    MeshLibTriLineMesh()
        : mesh(nullptr)
    {
        nodes.push_back(new MeshLib::Node(0, 0, 0));
        nodes.push_back(new MeshLib::Node(1, 0, 0));
        nodes.push_back(new MeshLib::Node(0, 1, 0));
        nodes.push_back(new MeshLib::Node(1, 1, 0));

        std::array<MeshLib::Node*, 3> t_nodes;
        t_nodes[0] = nodes[0];
        t_nodes[1] = nodes[1];
        t_nodes[2] = nodes[2];
        elements.push_back(new MeshLib::Tri(t_nodes));

        t_nodes[0] = nodes[3];
        t_nodes[1] = nodes[1];
        t_nodes[2] = nodes[2];
        elements.push_back(new MeshLib::Tri(t_nodes));

        std::array<MeshLib::Node*, 2> l_nodes;
        l_nodes[0] = nodes[1];
        l_nodes[1] = nodes[2];
        elements.push_back(new MeshLib::Line(l_nodes));

        mesh = new MeshLib::Mesh("M", nodes, elements);
    }

    ~MeshLibTriLineMesh()
    {
        /*std::remove_if(elements.begin(), elements.end(),
                [](MeshLib::Element* e) { delete e; return true; });
        std::remove_if(nodes.begin(), nodes.end(),
                [](MeshLib::Node* n) { delete n; return true; });
                */
        delete mesh;
    }

    MeshLib::Mesh const* mesh;
    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elements;

    public:
    // Helper functions to access elements and nodes.


    bool
    isConnectedToNode(std::size_t const n, std::size_t const e) const
    {
        return std::find(
            nodes[n]->getElements().cbegin(),
            nodes[n]->getElements().cend(),
            elements[e])
                != nodes[n]->getElements().cend();
    }
};

TEST_F(MeshLibTriLineMesh, Construction)
{
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_EQ(3u, mesh->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE, elements[0]->getGeomType());
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE, elements[1]->getGeomType());
    ASSERT_EQ(MeshLib::MeshElemType::LINE, elements[2]->getGeomType());
    ASSERT_EQ(4u, mesh->getNumberOfNodes());
}

TEST_F(MeshLibTriLineMesh, NodeToElementConnectivity)
{
    // Nodes 0 and 3 are connected only to triangles.
    EXPECT_EQ(1u, nodes[0]->getNumberOfElements());
    EXPECT_EQ(1u, nodes[3]->getNumberOfElements());
    EXPECT_EQ(elements[0], nodes[0]->getElement(0));
    EXPECT_EQ(elements[1], nodes[3]->getElement(0));

    // Nodes 1 and 2 are connected to all elements.
    EXPECT_EQ(3u, nodes[1]->getNumberOfElements());
    EXPECT_TRUE(isConnectedToNode(1, 0));
    EXPECT_TRUE(isConnectedToNode(1, 1));
    EXPECT_TRUE(isConnectedToNode(1, 2));

    EXPECT_EQ(3u, nodes[2]->getNumberOfElements());
    EXPECT_TRUE(isConnectedToNode(2, 0));
    EXPECT_TRUE(isConnectedToNode(2, 1));
    EXPECT_TRUE(isConnectedToNode(2, 2));
}

TEST_F(MeshLibTriLineMesh, ElementToElementConnectivity)
{
    // Triangles have other triangles as neighbors only.
    EXPECT_TRUE(elements[0]->hasNeighbor(elements[1]));
    EXPECT_TRUE(elements[1]->hasNeighbor(elements[0]));

    EXPECT_FALSE(elements[0]->hasNeighbor(elements[2]));
    EXPECT_FALSE(elements[1]->hasNeighbor(elements[2]));

    // Line has no neighbors
    EXPECT_FALSE(elements[2]->hasNeighbor(elements[0]));
    EXPECT_FALSE(elements[2]->hasNeighbor(elements[1]));
}
