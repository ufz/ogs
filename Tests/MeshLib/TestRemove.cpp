// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>

#include "MathLib/MathTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

TEST(MeshLib, RemoveNodes)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshToolsLib::MeshGenerator::generateLineMesh(1.0, 9)};

    std::vector<std::size_t> removed_node_ids;
    for (std::size_t i = 0; i < 5; i++)
    {
        removed_node_ids.push_back(i);
    }

    auto new_mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshToolsLib::removeNodes(*mesh, removed_node_ids, "")};

    ASSERT_EQ(5u, new_mesh->getNumberOfNodes());
    ASSERT_EQ(5u, new_mesh->computeNumberOfBaseNodes());
    ASSERT_EQ(4u, new_mesh->getNumberOfElements());
    for (std::size_t i = 0; i < new_mesh->getNumberOfNodes(); i++)
    {
        ASSERT_TRUE(*mesh->getNode(5 + i) == *new_mesh->getNode(i));
    }
}

TEST(MeshLib, RemoveElements)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshToolsLib::MeshGenerator::generateLineMesh(1.0, 9)};

    std::vector<std::size_t> removed_ele_ids;
    for (std::size_t i = 0; i < 5; i++)
    {
        removed_ele_ids.push_back(i);
    }

    auto new_mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshToolsLib::removeElements(*mesh, removed_ele_ids, "")};

    ASSERT_EQ(5u, new_mesh->getNumberOfNodes());
    ASSERT_EQ(5u, new_mesh->computeNumberOfBaseNodes());
    ASSERT_EQ(4u, new_mesh->getNumberOfElements());
    for (std::size_t i = 0; i < new_mesh->getNumberOfNodes(); i++)
    {
        ASSERT_TRUE(*mesh->getNode(5 + i) == *new_mesh->getNode(i));
    }
}
