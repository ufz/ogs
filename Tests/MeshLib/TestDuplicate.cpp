/**
 * @file TestDuplicate.cpp
 * @author Karsten Rink
 * @date 2013-03-25
 * @brief Tests for Duplicate functions
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include "MathLib/MathTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshQuality/MeshValidation.h"
#include "MeshLib/Node.h"

TEST(MeshLib, Duplicate)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshLib::MeshGenerator::generateRegularQuadMesh(10, 5, 1)};

    std::vector<MeshLib::Node*> new_nodes (MeshLib::copyNodeVector(mesh->getNodes()));
    std::vector<MeshLib::Element*> new_elements (MeshLib::copyElementVector(mesh->getElements(), new_nodes));

    MeshLib::Mesh new_mesh ("new", new_nodes, new_elements);

    ASSERT_EQ (mesh->getNumberOfElements(), new_mesh.getNumberOfElements());
    ASSERT_EQ (mesh->getNumberOfNodes(), new_mesh.getNumberOfNodes());

    std::vector<std::size_t> del_idx(1,1);
    std::unique_ptr<MeshLib::Mesh> mesh2(MeshLib::removeNodes(*mesh, del_idx, "mesh2"));

    ASSERT_EQ (mesh2->getNumberOfElements(), new_mesh.getNumberOfElements()-2);
    ASSERT_EQ (mesh2->getNumberOfNodes(), new_mesh.getNumberOfNodes()-2);

    ASSERT_DOUBLE_EQ (4.0, MathLib::sqrDist(*mesh2->getNode(0), *new_mesh.getNode(0)));
    ASSERT_DOUBLE_EQ (0.0, MathLib::sqrDist(*mesh2->getNode(0), *new_mesh.getNode(2)));

    ASSERT_DOUBLE_EQ (4.0, MathLib::sqrDist(*mesh2->getElement(0)->getNode(0), *new_mesh.getElement(0)->getNode(0)));
    ASSERT_DOUBLE_EQ (0.0, MathLib::sqrDist(*mesh2->getElement(0)->getNode(0), *new_mesh.getElement(2)->getNode(0)));
}
