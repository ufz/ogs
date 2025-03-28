/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <memory>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/FlipElements.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

TEST(MeshLib, FlipLineMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshToolsLib::MeshGenerator::generateLineMesh(1.0, 5));
    std::unique_ptr<MeshLib::Mesh> result(
        MeshToolsLib::createFlippedMesh(*mesh));

    ASSERT_EQ(mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(mesh->getNumberOfElements(), result->getNumberOfElements());
    for (std::size_t i = 0; i < result->getNumberOfElements(); ++i)
    {
        ASSERT_EQ(mesh->getElement(i)->getNode(0)->getID(),
                  result->getElement(i)->getNode(1)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(1)->getID(),
                  result->getElement(i)->getNode(0)->getID());
    }
}

TEST(MeshLib, FlipTriMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshToolsLib::MeshGenerator::generateRegularTriMesh(5, 5));
    std::unique_ptr<MeshLib::Mesh> result(
        MeshToolsLib::createFlippedMesh(*mesh));

    ASSERT_EQ(mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(mesh->getNumberOfElements(), result->getNumberOfElements());
    for (std::size_t i = 0; i < result->getNumberOfElements(); ++i)
    {
        ASSERT_EQ(mesh->getElement(i)->getNode(0)->getID(),
                  result->getElement(i)->getNode(1)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(1)->getID(),
                  result->getElement(i)->getNode(0)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(2)->getID(),
                  result->getElement(i)->getNode(2)->getID());
        ASSERT_EQ(
            1.0,
            MeshLib::FaceRule::getSurfaceNormal(*result->getElement(i))[2]);
    }
}

TEST(MeshLib, FlipQuadMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshToolsLib::MeshGenerator::generateRegularQuadMesh(5, 5));
    std::unique_ptr<MeshLib::Mesh> result(
        MeshToolsLib::createFlippedMesh(*mesh));

    ASSERT_EQ(mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(mesh->getNumberOfElements(), result->getNumberOfElements());
    for (std::size_t i = 0; i < result->getNumberOfElements(); ++i)
    {
        ASSERT_EQ(mesh->getElement(i)->getNode(0)->getID(),
                  result->getElement(i)->getNode(1)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(1)->getID(),
                  result->getElement(i)->getNode(0)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(2)->getID(),
                  result->getElement(i)->getNode(3)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(3)->getID(),
                  result->getElement(i)->getNode(2)->getID());
        ASSERT_EQ(
            1.0,
            MeshLib::FaceRule::getSurfaceNormal(*result->getElement(i))[2]);
    }
}

TEST(MeshLib, FlipHexMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshToolsLib::MeshGenerator::generateRegularHexMesh(2, 2));
    std::unique_ptr<MeshLib::Mesh> result(
        MeshToolsLib::createFlippedMesh(*mesh));

    ASSERT_EQ(nullptr, result);
    std::vector<MeshLib::Node*> nodes;
    for (std::size_t i = 0; i < mesh->getNumberOfNodes(); ++i)
    {
        nodes.push_back(new MeshLib::Node(*mesh->getNode(i)));
    }
    std::unique_ptr<MeshLib::Element> elem(
        MeshToolsLib::createFlippedElement(*mesh->getElement(0), nodes));
    ASSERT_EQ(nullptr, elem);
    for (MeshLib::Node* n : nodes)
    {
        delete n;
    }
}

TEST(MeshLib, DoubleFlipQuadMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshToolsLib::MeshGenerator::generateRegularQuadMesh(5, 5));
    std::unique_ptr<MeshLib::Mesh> result(
        MeshToolsLib::createFlippedMesh(*mesh));
    std::unique_ptr<MeshLib::Mesh> result2(
        MeshToolsLib::createFlippedMesh(*result));

    ASSERT_EQ(mesh->getNumberOfNodes(), result2->getNumberOfNodes());
    ASSERT_EQ(mesh->getNumberOfElements(), result2->getNumberOfElements());
    for (std::size_t i = 0; i < result2->getNumberOfElements(); ++i)
    {
        for (std::size_t j = 0; j < 4; ++j)
        {
            ASSERT_EQ(mesh->getElement(i)->getNode(j)->getID(),
                      result2->getElement(i)->getNode(j)->getID());
        }
    }
}
