/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Tet.h"

using namespace MeshLib;

TEST(MeshLib, ElementConstantsQuad4)
{
    ASSERT_EQ(2u, Quad::dimension);
    ASSERT_EQ(4u, Quad::n_all_nodes);
    ASSERT_EQ(4u, Quad::n_base_nodes);

    std::array<Node*, 4> nodes;
    nodes[0] = new Node(0.0, 0.0, 0.0);
    nodes[1] = new Node(0.0, 1.0, 0.0);
    nodes[2] = new Node(1.0, 1.0, 0.0);
    nodes[3] = new Node(1.0, 0.0, 0.0);
    Quad quad(nodes);

    ASSERT_EQ(2u, quad.getDimension());
    ASSERT_EQ(4u, quad.getNumberOfNodes());
    ASSERT_EQ(4u, quad.getNumberOfBaseNodes());

    for (auto n : nodes)
        delete n;
}

TEST(MeshLib, ElementConstantsQuad8)
{
    ASSERT_EQ(2u, Quad8::dimension);
    ASSERT_EQ(8u, Quad8::n_all_nodes);
    ASSERT_EQ(4u, Quad8::n_base_nodes);

    std::array<Node*, 8> nodes;
    nodes[0] = new Node(0.0, 0.0, 0.0);
    nodes[1] = new Node(0.0, 1.0, 0.0);
    nodes[2] = new Node(1.0, 1.0, 0.0);
    nodes[3] = new Node(1.0, 0.0, 0.0);

    nodes[4] = new Node(0.5, 0.0, 0.0);
    nodes[5] = new Node(1.0, 0.5, 0.0);
    nodes[6] = new Node(0.5, 1.0, 0.0);
    nodes[7] = new Node(0.5, 0.0, 0.0);

    Quad8 quad8(nodes);

    ASSERT_EQ(2u, quad8.getDimension());
    ASSERT_EQ(8u, quad8.getNumberOfNodes());
    ASSERT_EQ(4u, quad8.getNumberOfBaseNodes());

    for (auto n : nodes)
        delete n;
}

TEST(MeshLib, ElementConstantsQuad9)
{
    ASSERT_EQ(2u, Quad9::dimension);
    ASSERT_EQ(9u, Quad9::n_all_nodes);
    ASSERT_EQ(4u, Quad9::n_base_nodes);

    std::array<Node*, 9> nodes;
    nodes[0] = new Node(0.0, 0.0, 0.0);
    nodes[1] = new Node(0.0, 1.0, 0.0);
    nodes[2] = new Node(1.0, 1.0, 0.0);
    nodes[3] = new Node(1.0, 0.0, 0.0);

    nodes[4] = new Node(0.5, 0.0, 0.0);
    nodes[5] = new Node(1.0, 0.5, 0.0);
    nodes[6] = new Node(0.5, 1.0, 0.0);
    nodes[7] = new Node(0.5, 0.0, 0.0);

    nodes[8] = new Node(0.5, 0.5, 0.0);
    Quad9 quad9(nodes);

    ASSERT_EQ(2u, quad9.getDimension());
    ASSERT_EQ(9u, quad9.getNumberOfNodes());
    ASSERT_EQ(4u, quad9.getNumberOfBaseNodes());

    for (auto n : nodes)
        delete n;
}

TEST(MeshLib, ElementConstantsHex8)
{
    ASSERT_EQ(3u, Hex::dimension);
    ASSERT_EQ(8u, Hex::n_all_nodes);
    ASSERT_EQ(8u, Hex::n_base_nodes);

    std::array<Node*, 8> nodes;
    nodes[0] = new Node(0.0, 0.0, 0.0);
    nodes[1] = new Node(0.0, 1.0, 0.0);
    nodes[2] = new Node(1.0, 1.0, 0.0);
    nodes[3] = new Node(1.0, 0.0, 0.0);
    nodes[4] = new Node(0.0, 0.0, 1.0);
    nodes[5] = new Node(0.0, 1.0, 1.0);
    nodes[6] = new Node(1.0, 1.0, 1.0);
    nodes[7] = new Node(1.0, 0.0, 1.0);
    Hex ele(nodes);

    ASSERT_EQ(3u, ele.getDimension());
    ASSERT_EQ(8u, ele.getNumberOfNodes());
    ASSERT_EQ(8u, ele.getNumberOfBaseNodes());

    for (auto n : nodes)
        delete n;
}

TEST(MeshLib, ElementConstantsHex20)
{
    ASSERT_EQ( 3u, Hex20::dimension);
    ASSERT_EQ(20u, Hex20::n_all_nodes);
    ASSERT_EQ( 8u, Hex20::n_base_nodes);

    std::array<Node*, 20> nodes;
    nodes[0] = new Node(0.0, 0.0, 0.0);
    nodes[1] = new Node(1.0, 0.0, 0.0);
    nodes[2] = new Node(1.0, 1.0, 0.0);
    nodes[3] = new Node(0.0, 1.0, 0.0);
    nodes[4] = new Node(0.0, 0.0, 1.0);
    nodes[5] = new Node(1.0, 0.0, 1.0);
    nodes[6] = new Node(1.0, 1.0, 1.0);
    nodes[7] = new Node(0.0, 1.0, 1.0);
    nodes[8] = new Node(0.5, 0.0, 0.0);
    nodes[9] = new Node(1.0, 0.5, 0.0);
    nodes[10] = new Node(0.5, 1.0, 0.0);
    nodes[11] = new Node(0.0, 0.5, 0.0);
    nodes[12] = new Node(0.5, 0.0, 1.0);
    nodes[13] = new Node(1.0, 0.5, 1.0);
    nodes[14] = new Node(0.5, 1.0, 1.0);
    nodes[15] = new Node(0.0, 0.5, 1.0);
    nodes[16] = new Node(0.0, 0.0, 0.5);
    nodes[17] = new Node(1.0, 0.0, 0.5);
    nodes[18] = new Node(1.0, 1.0, 0.5);
    nodes[19] = new Node(0.0, 1.0, 0.5);
    Hex20 ele(nodes);

    ASSERT_EQ( 3u, ele.getDimension());
    ASSERT_EQ(20u, ele.getNumberOfNodes());
    ASSERT_EQ( 8u, ele.getNumberOfBaseNodes());

    for (auto n : nodes)
        delete n;
}

TEST(MeshLib, ElementConstantsTet4)
{
    ASSERT_EQ(3u, Tet::dimension);
    ASSERT_EQ(4u, Tet::n_all_nodes);
    ASSERT_EQ(4u, Tet::n_base_nodes);

    std::array<Node*, 4> nodes;
    nodes[0] = new Node(0.0, 0.0, 0.0);
    nodes[1] = new Node(1.0, 0.0, 0.0);
    nodes[2] = new Node(0.0, 1.0, 0.0);
    nodes[3] = new Node(0.0, 0.0, 1.0);
    Tet ele(nodes);

    ASSERT_EQ(3u, ele.getDimension());
    ASSERT_EQ(4u, ele.getNumberOfNodes());
    ASSERT_EQ(4u, ele.getNumberOfBaseNodes());

    for (auto n : nodes)
        delete n;
}

TEST(MeshLib, ElementConstantsTet10)
{
    ASSERT_EQ( 3u, Tet10::dimension);
    ASSERT_EQ(10u, Tet10::n_all_nodes);
    ASSERT_EQ( 4u, Tet10::n_base_nodes);

    std::array<Node*, 10> nodes;
    nodes[0] = new Node(0.0, 0.0, 0.0);
    nodes[1] = new Node(1.0, 0.0, 0.0);
    nodes[2] = new Node(0.0, 1.0, 0.0);
    nodes[3] = new Node(0.0, 0.0, 1.0);

    nodes[4] = new Node(0.5, 0.0, 0.0);
    nodes[5] = new Node(0.5, 0.5, 0.0);
    nodes[6] = new Node(0.0, 0.5, 0.0);
    nodes[7] = new Node(0.0, 0.0, 0.5);
    nodes[8] = new Node(0.5, 0.0, 0.5);
    nodes[9] = new Node(0.0, 0.5, 0.5);
    Tet10 ele(nodes);

    ASSERT_EQ( 3u, ele.getDimension());
    ASSERT_EQ(10u, ele.getNumberOfNodes());
    ASSERT_EQ( 4u, ele.getNumberOfBaseNodes());

    for (auto n : nodes)
        delete n;
}

