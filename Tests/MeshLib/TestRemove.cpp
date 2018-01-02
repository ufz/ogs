/**
 * @file TestDuplicate.cpp
 * @author Karsten Rink
 * @date 2013-03-25
 * @brief Tests for Duplicate functions
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include "MathLib/MathTools.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshQuality/MeshValidation.h"

TEST(MeshLib, RemoveNodes)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshLib::MeshGenerator::generateLineMesh(1.0, 9)};

    std::vector<std::size_t> removed_node_ids;
    for (std::size_t i=0; i<5; i++)
        removed_node_ids.push_back(i);

    auto new_mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshLib::removeNodes(*mesh, removed_node_ids, "")};

    ASSERT_EQ(5u, new_mesh->getNumberOfNodes());
    ASSERT_EQ(5u, new_mesh->getNumberOfBaseNodes());
    ASSERT_EQ(4u, new_mesh->getNumberOfElements());
    for (std::size_t i=0; i<new_mesh->getNumberOfNodes(); i++)
        ASSERT_TRUE(*mesh->getNode(5+i) == *new_mesh->getNode(i));
}

TEST(MeshLib, RemoveElements)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshLib::MeshGenerator::generateLineMesh(1.0, 9)};

    std::vector<std::size_t> removed_ele_ids;
    for (std::size_t i=0; i<5; i++)
        removed_ele_ids.push_back(i);

    auto new_mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshLib::removeElements(*mesh, removed_ele_ids, "")};

    ASSERT_EQ(5u, new_mesh->getNumberOfNodes());
    ASSERT_EQ(5u, new_mesh->getNumberOfBaseNodes());
    ASSERT_EQ(4u, new_mesh->getNumberOfElements());
    for (std::size_t i=0; i<new_mesh->getNumberOfNodes(); i++)
        ASSERT_TRUE(*mesh->getNode(5+i) == *new_mesh->getNode(i));
}
