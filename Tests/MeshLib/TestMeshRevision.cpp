// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/MeshRevision.h"

TEST(MeshEditing, Tri)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0.1));

    std::array<MeshLib::Node*, 3> nodes_array = {
        {nodes[0], nodes[1], nodes[2]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Tri(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(MeshLib::MeshElemType::LINE,
              result->getElement(0)->getGeomType());
    ASSERT_EQ(1, result->getElement(0)->getContent());
    ASSERT_EQ(2u, result->getNumberOfNodes());
    delete result;

    result = rev.simplifyMesh("new_mesh", 0.0999);
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE,
              result->getElement(0)->getGeomType());
    ASSERT_EQ(0.05, result->getElement(0)->getContent());
    ASSERT_EQ(3u, result->getNumberOfNodes());
    delete result;
}

TEST(MeshEditing, NodePropertyAfterCollapse)
{
    // Triangle whose nodes 1 and 2 are within eps and get collapsed, so the
    // three input nodes are reduced to two output nodes. The surviving node is
    // always the one with the lower index (node 2 collapses into node 1).
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0.1));

    std::array<MeshLib::Node*, 3> const nodes_array = {
        {nodes[0], nodes[1], nodes[2]}};
    std::vector<MeshLib::Element*> elements;
    elements.push_back(new MeshLib::Tri(nodes_array));
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    // One property per value type supported by the revision, to cover all of
    // them for node properties.
    auto* const double_prop =
        mesh.getProperties().createNewPropertyVector<double>(
            "NodeValuesDouble", MeshLib::MeshItemType::Node, 1);
    ASSERT_NE(nullptr, double_prop);
    double_prop->assign(std::vector<double>{100.0, 200.0, 300.0});

    auto* const float_prop =
        mesh.getProperties().createNewPropertyVector<float>(
            "NodeValuesFloat", MeshLib::MeshItemType::Node, 1);
    ASSERT_NE(nullptr, float_prop);
    float_prop->assign(std::vector<float>{1.5f, 2.5f, 3.5f});

    auto* const int_prop = mesh.getProperties().createNewPropertyVector<int>(
        "NodeValuesInt", MeshLib::MeshItemType::Node, 1);
    ASSERT_NE(nullptr, int_prop);
    int_prop->assign(std::vector<int>{10, 20, 30});

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* const result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfNodes());

    auto const* new_double_prop =
        result->getProperties().getPropertyVector<double>("NodeValuesDouble");
    ASSERT_NE(nullptr, new_double_prop);
    // The property must be compacted to match the reduced node array: one value
    // per surviving node, in the new node order (no gaps, no stale size).
    ASSERT_EQ(result->getNumberOfNodes(), new_double_prop->getNumberOfTuples());
    ASSERT_EQ(100.0, (*new_double_prop)[0]);  // original node 0
    ASSERT_EQ(200.0,
              (*new_double_prop)[1]);  // original node 1 (node 2 collapsed
                                       // here)

    auto const* new_float_prop =
        result->getProperties().getPropertyVector<float>("NodeValuesFloat");
    ASSERT_NE(nullptr, new_float_prop);
    ASSERT_EQ(result->getNumberOfNodes(), new_float_prop->getNumberOfTuples());
    ASSERT_EQ(1.5f, (*new_float_prop)[0]);
    ASSERT_EQ(2.5f, (*new_float_prop)[1]);

    auto const* new_int_prop =
        result->getProperties().getPropertyVector<int>("NodeValuesInt");
    ASSERT_NE(nullptr, new_int_prop);
    ASSERT_EQ(result->getNumberOfNodes(), new_int_prop->getNumberOfTuples());
    ASSERT_EQ(10, (*new_int_prop)[0]);
    ASSERT_EQ(20, (*new_int_prop)[1]);

    delete result;
}

TEST(MeshEditing, CellPropertyAfterSubdivision)
{
    // Non-planar quad which is subdivided into two triangles, i.e. the revised
    // mesh has more elements than the input mesh. Every new element takes the
    // value of the element it originates from.
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0.1));
    nodes.push_back(new MeshLib::Node(1, 0, 0));

    std::array<MeshLib::Node*, 4> const nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3]}};
    std::vector<MeshLib::Element*> elements;
    elements.push_back(new MeshLib::Quad(nodes_array));
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    // One property per value type supported by the revision, to cover all of
    // them for cell properties.
    auto* const int_prop = mesh.getProperties().createNewPropertyVector<int>(
        "CellValuesInt", MeshLib::MeshItemType::Cell, 1);
    ASSERT_NE(nullptr, int_prop);
    int_prop->assign(std::vector<int>{42});

    auto* const float_prop =
        mesh.getProperties().createNewPropertyVector<float>(
            "CellValuesFloat", MeshLib::MeshItemType::Cell, 1);
    ASSERT_NE(nullptr, float_prop);
    float_prop->assign(std::vector<float>{4.25f});

    auto* const double_prop =
        mesh.getProperties().createNewPropertyVector<double>(
            "CellValuesDouble", MeshLib::MeshItemType::Cell, 1);
    ASSERT_NE(nullptr, double_prop);
    double_prop->assign(std::vector<double>{13.75});

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* const result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());

    auto const* new_int_prop =
        result->getProperties().getPropertyVector<int>("CellValuesInt");
    ASSERT_NE(nullptr, new_int_prop);
    // One value per element of the revised mesh, not per element of the input
    // mesh.
    ASSERT_EQ(result->getNumberOfElements(), new_int_prop->getNumberOfTuples());
    ASSERT_EQ(42, (*new_int_prop)[0]);
    ASSERT_EQ(42, (*new_int_prop)[1]);

    auto const* new_float_prop =
        result->getProperties().getPropertyVector<float>("CellValuesFloat");
    ASSERT_NE(nullptr, new_float_prop);
    ASSERT_EQ(result->getNumberOfElements(),
              new_float_prop->getNumberOfTuples());
    ASSERT_EQ(4.25f, (*new_float_prop)[0]);
    ASSERT_EQ(4.25f, (*new_float_prop)[1]);

    auto const* new_double_prop =
        result->getProperties().getPropertyVector<double>("CellValuesDouble");
    ASSERT_NE(nullptr, new_double_prop);
    ASSERT_EQ(result->getNumberOfElements(),
              new_double_prop->getNumberOfTuples());
    ASSERT_EQ(13.75, (*new_double_prop)[0]);
    ASSERT_EQ(13.75, (*new_double_prop)[1]);

    delete result;
}

TEST(MeshEditing, NonPlanarQuad)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0.1));
    nodes.push_back(new MeshLib::Node(1, 0, 0));

    std::array<MeshLib::Node*, 4> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Quad(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);
    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE,
              result->getElement(1)->getGeomType());

    delete result;
}

TEST(MeshEditing, Quad2Line)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0.1));
    nodes.push_back(new MeshLib::Node(1, 0, 0.1));

    std::array<MeshLib::Node*, 4> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Quad(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(MeshLib::MeshElemType::LINE,
              result->getElement(0)->getGeomType());
    ASSERT_NEAR(1.414213562373095, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_EQ(2u, result->getNumberOfNodes());

    delete result;
}

TEST(MeshEditing, Quad2Tri)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0.1));
    nodes.push_back(new MeshLib::Node(1, 1, 0.1));

    std::array<MeshLib::Node*, 4> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Quad(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE,
              result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.5049752469181039, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_EQ(3u, result->getNumberOfNodes());

    delete result;
}

TEST(MeshEditing, NonPlanarHex)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0, 0, -0.5));
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 1));
    nodes.push_back(new MeshLib::Node(1, 0, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(0, 1, 1));

    std::array<MeshLib::Node*, 8> nodes_array = {{nodes[0], nodes[1], nodes[2],
                                                  nodes[3], nodes[4], nodes[5],
                                                  nodes[6], nodes[7]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Hex(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(6u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
              result->getElement(4)->getGeomType());
    ASSERT_NEAR(0.25, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(5)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Hex2PyramidPrism)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 1));
    nodes.push_back(new MeshLib::Node(1, 0, .1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(0, 1, 1));

    std::array<MeshLib::Node*, 8> nodes_array = {{nodes[0], nodes[1], nodes[2],
                                                  nodes[3], nodes[4], nodes[5],
                                                  nodes[6], nodes[7]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Hex(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::PYRAMID,
              result->getElement(0)->getGeomType());
    ASSERT_EQ(MeshLib::MeshElemType::PRISM,
              result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.3333333333333333, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.5, result->getElement(1)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Hex2FourTets)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 1));
    nodes.push_back(new MeshLib::Node(1, 0, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(0, 1, 1));

    std::array<MeshLib::Node*, 8> nodes_array = {{nodes[0], nodes[1], nodes[2],
                                                  nodes[3], nodes[4], nodes[5],
                                                  nodes[6], nodes[7]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Hex(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(4u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
              result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(1)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(2)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(3)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Hex2TwoTets)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));

    std::array<MeshLib::Node*, 8> nodes_array = {{nodes[0], nodes[1], nodes[2],
                                                  nodes[3], nodes[4], nodes[5],
                                                  nodes[6], nodes[7]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Hex(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
              result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(1)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, NonPlanarPyramid)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0, -.5));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 1));

    std::array<MeshLib::Node*, 5> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Pyramid(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
              result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.25, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.1666666666666667, result->getElement(1)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Pyramid2Tet)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 1));

    std::array<MeshLib::Node*, 5> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Pyramid(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
              result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.16666666666666666, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Pyramid2Quad)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 0.1));

    std::array<MeshLib::Node*, 5> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Pyramid(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::QUAD,
              result->getElement(0)->getGeomType());
    ASSERT_NEAR(1, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Pyramid2Tri)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0.1, 0));
    nodes.push_back(new MeshLib::Node(0, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 0.1));

    std::array<MeshLib::Node*, 5> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Pyramid(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(3u, result->getNumberOfNodes());
    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE,
              result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.5, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, NonPlanarPrism)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(0, -0.5, 2));

    std::array<MeshLib::Node*, 6> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(3u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
              result->getElement(2)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Prism2TwoTets)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0.9, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(0, 0, 1));

    std::array<MeshLib::Node*, 6> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(5u, result->getNumberOfNodes());
    ASSERT_EQ(2u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
              result->getElement(1)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.15, result->getElement(1)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Prism2Quad)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0.9, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0.9, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(0, 0, 1));

    std::array<MeshLib::Node*, 6> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::QUAD,
              result->getElement(0)->getGeomType());
    ASSERT_NEAR(1.345362404707371, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Prism2Tet)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(0.9, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 0.9, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(0, 0, 1));

    std::array<MeshLib::Node*, 6> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
              result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.1666666666666667, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

TEST(MeshEditing, Prism2Tri)
{
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(1, 0, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 0));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(1, 1, 1));
    nodes.push_back(new MeshLib::Node(0.9, 0.9, 1));

    std::array<MeshLib::Node*, 6> nodes_array = {
        {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]}};
    std::vector<MeshLib::Element*> elements;
    MeshLib::Element* elem(new MeshLib::Prism(nodes_array));
    elements.push_back(elem);
    MeshLib::Mesh mesh("testmesh", nodes, elements);

    MeshToolsLib::MeshRevision rev(mesh);
    MeshLib::Mesh* result = rev.simplifyMesh("new_mesh", 0.2);

    ASSERT_EQ(1u, result->getNumberOfElements());
    ASSERT_EQ(MeshLib::MeshElemType::TRIANGLE,
              result->getElement(0)->getGeomType());
    ASSERT_NEAR(0.5, result->getElement(0)->getContent(),
                std::numeric_limits<double>::epsilon());

    delete result;
}

namespace MeshToolsLib
{
extern unsigned lutPrismThirdNode(unsigned, unsigned);
}

TEST(MeshEditing, PrismLutThirdNode)
{
    constexpr auto max = std::numeric_limits<unsigned>::max();
    // clang-format off
    // Expected results in a matrix 6x6 to cover all possible combinations
    constexpr std::array<std::array<unsigned, 6>, 6> expected =
        //   0    1    2    3    4    5
        {{{{max, 2  , 1  , max, max, max}},     // 0
          {{2  , max, 0  , max, max, max}},     // 1
          {{1  , 0  , max, max, max, max}},     // 2
          {{max, max, max, max, 5  , 4  }},     // 3
          {{max, max, max, 5  , max, 3  }},     // 4
          {{max, max, max, 4  , 3  , max}}}};   // 5
    // clang-format on

    for (unsigned i = 0; i < 6; ++i)
    {
        for (unsigned j = 0; j < 6; ++j)
        {
            EXPECT_EQ(MeshToolsLib::lutPrismThirdNode(i, j), expected[i][j])
                << "for (i,j) = " << i << ", " << j;
        }
    }

    // Test for outside values
    EXPECT_EQ(MeshToolsLib::lutPrismThirdNode(-1, -1), max);
    EXPECT_EQ(MeshToolsLib::lutPrismThirdNode(0, -1), max);
    EXPECT_EQ(MeshToolsLib::lutPrismThirdNode(-1, 0), max);
    EXPECT_EQ(MeshToolsLib::lutPrismThirdNode(5, 6), max);
    EXPECT_EQ(MeshToolsLib::lutPrismThirdNode(6, 5), max);
    EXPECT_EQ(MeshToolsLib::lutPrismThirdNode(6, 6), max);
}
