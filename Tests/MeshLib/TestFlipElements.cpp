/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshEditing/FlipElements.h"

TEST(MeshLib, FlipLineMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::MeshGenerator::generateLineMesh(1.0, 5));
    std::unique_ptr<MeshLib::Mesh> result (MeshLib::createFlippedMesh(*mesh));

    ASSERT_EQ(mesh->getNNodes(), result->getNNodes());
    ASSERT_EQ(mesh->getNElements(), result->getNElements());
    for (std::size_t i=0; i<result->getNElements(); ++i)
    {
        ASSERT_EQ(mesh->getElement(i)->getNode(0)->getID(),
                  result->getElement(i)->getNode(1)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(1)->getID(),
                  result->getElement(i)->getNode(0)->getID());
    }
}

TEST(MeshLib, FlipTriMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::MeshGenerator::generateRegularTriMesh(5, 5));
    std::unique_ptr<MeshLib::Mesh> result (MeshLib::createFlippedMesh(*mesh));

    ASSERT_EQ(mesh->getNNodes(), result->getNNodes());
    ASSERT_EQ(mesh->getNElements(), result->getNElements());
    for (std::size_t i=0; i<result->getNElements(); ++i)
    {
        ASSERT_EQ(mesh->getElement(i)->getNode(0)->getID(),
                  result->getElement(i)->getNode(1)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(1)->getID(),
                  result->getElement(i)->getNode(0)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(2)->getID(),
                  result->getElement(i)->getNode(2)->getID());
        ASSERT_EQ(1.0, MeshLib::FaceRule::getSurfaceNormal(result->getElement(i))[2]);
    }
}

TEST(MeshLib, FlipQuadMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::MeshGenerator::generateRegularQuadMesh(5, 5));
    std::unique_ptr<MeshLib::Mesh> result (MeshLib::createFlippedMesh(*mesh));

    ASSERT_EQ(mesh->getNNodes(), result->getNNodes());
    ASSERT_EQ(mesh->getNElements(), result->getNElements());
    for (std::size_t i=0; i<result->getNElements(); ++i)
    {
        ASSERT_EQ(mesh->getElement(i)->getNode(0)->getID(),
                  result->getElement(i)->getNode(1)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(1)->getID(),
                  result->getElement(i)->getNode(0)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(2)->getID(),
                  result->getElement(i)->getNode(3)->getID());
        ASSERT_EQ(mesh->getElement(i)->getNode(3)->getID(),
                  result->getElement(i)->getNode(2)->getID());
        ASSERT_EQ(1.0, MeshLib::FaceRule::getSurfaceNormal(result->getElement(i))[2]);
    }
}

TEST(MeshLib, FlipHexMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::MeshGenerator::generateRegularHexMesh(2, 2));
    std::unique_ptr<MeshLib::Mesh> result (MeshLib::createFlippedMesh(*mesh));

    ASSERT_EQ(nullptr, result);
    std::vector<MeshLib::Node*> nodes;
    for (std::size_t i=0; i<mesh->getNNodes(); ++i)
        nodes.push_back(new MeshLib::Node(*mesh->getNode(i)));
    std::unique_ptr<MeshLib::Element> elem (MeshLib::createFlippedElement(*mesh->getElement(0), nodes));
    ASSERT_EQ(nullptr, elem);
    for (MeshLib::Node* n : nodes)
        delete n;
}

TEST(MeshLib, DoubleFlipQuadMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::MeshGenerator::generateRegularQuadMesh(5, 5));
    std::unique_ptr<MeshLib::Mesh> result (MeshLib::createFlippedMesh(*mesh));
    std::unique_ptr<MeshLib::Mesh> result2 (MeshLib::createFlippedMesh(*result));

    ASSERT_EQ(mesh->getNNodes(), result2->getNNodes());
    ASSERT_EQ(mesh->getNElements(), result2->getNElements());
    for (std::size_t i=0; i<result2->getNElements(); ++i)
        for (std::size_t j=0; j<4; ++j)
            ASSERT_EQ(mesh->getElement(i)->getNode(j)->getID(),
                      result2->getElement(i)->getNode(j)->getID());
}

