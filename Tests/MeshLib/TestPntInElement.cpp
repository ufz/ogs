/**
 * @file TestPntInElement.cpp
 * @author Karsten Rink
 * @date 2014-09-23
 * @brief Tests for check if a point is located inside an element
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "Mesh.h"
#include "Node.h"
#include "MeshEditing/MeshRevision.h"
#include "Elements/Element.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"


std::vector<MeshLib::Node*> createNodes()
{
	std::vector<MeshLib::Node*> nodes;
	nodes.push_back(new MeshLib::Node(0,0,0));
	nodes.push_back(new MeshLib::Node(1,0,0));
	nodes.push_back(new MeshLib::Node(1,1,0));
    nodes.push_back(new MeshLib::Node(0,1,0));
	nodes.push_back(new MeshLib::Node(0,0,1));
	nodes.push_back(new MeshLib::Node(1,0,1));
	nodes.push_back(new MeshLib::Node(1,1,1));
    nodes.push_back(new MeshLib::Node(0,1,1));
    return nodes;
}

TEST(IsPntInElement, Line)
{
    GeoLib::Point pnt;
    std::vector<MeshLib::Node*> nodes (createNodes());
    std::array<MeshLib::Node*, 2> line_nodes = { nodes[0], nodes[4] };
    MeshLib::Line line(line_nodes);
    pnt = GeoLib::Point(0,0,0.7);
    ASSERT_EQ(true, line.isPntInElement(pnt));
    pnt = GeoLib::Point(0,0.1,0.7);
    ASSERT_EQ(false, line.isPntInElement(pnt));
}
TEST(IsPntInElement, Tri)
{
    GeoLib::Point pnt;
    std::vector<MeshLib::Node*> nodes (createNodes());
    std::array<MeshLib::Node*, 3> tri_nodes = { nodes[0], nodes[1], nodes[4] };
    MeshLib::Tri tri(tri_nodes);

    pnt = GeoLib::Point(0.1,0,.1);
    ASSERT_EQ(true, tri.isPntInElement(pnt));
    pnt = GeoLib::Point(0.9,0,.7);
    ASSERT_EQ(false, tri.isPntInElement(pnt));

}

TEST(IsPntInElement, Tet)
{
    GeoLib::Point pnt;
    std::vector<MeshLib::Node*> nodes (createNodes());
    std::array<MeshLib::Node*, 4> tet_nodes = { nodes[0], nodes[1], nodes[2], nodes[4] };
    MeshLib::Tet tet(tet_nodes);

    pnt = GeoLib::Point(0.5,0.3,0.1);
    ASSERT_EQ(true, tet.isPntInElement(pnt));
    pnt = GeoLib::Point(0.5,0.6,0.1);
    ASSERT_EQ(false, tet.isPntInElement(pnt));
}

TEST(IsPntInElement, Hex)
{
    GeoLib::Point pnt;
    std::vector<MeshLib::Node*> nodes (createNodes());
    std::array<MeshLib::Node*, 8> hex_nodes;
    std::copy(nodes.begin(), nodes.end(), hex_nodes.begin());
    MeshLib::Hex hex(hex_nodes);

    pnt = GeoLib::Point(0.99,0.99,0.99);
    ASSERT_EQ(true, hex.isPntInElement(pnt));
    pnt = GeoLib::Point(0.99,0,0);
    ASSERT_EQ(true, hex.isPntInElement(pnt));
    pnt = GeoLib::Point(0.0, 0.0, 0.0);
    ASSERT_EQ(true, hex.isPntInElement(pnt));
    pnt = GeoLib::Point(1.01,0.99,0.99);
    ASSERT_EQ(false, hex.isPntInElement(pnt));
}