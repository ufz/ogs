/**
 * @file TestProjectMeshOnPlane.cpp
 * @author Karsten Rink
 * @date 2015-04-16
 * @brief Tests for projectMeshOnPlane
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/projectMeshOntoPlane.h"
#include "MeshLib/Node.h"

class ProjectionTest : public ::testing::Test
{
public:
    ProjectionTest()
    : _mesh(nullptr)
    {
        std::size_t const n_nodes (100);
        std::vector<MeshLib::Node*> nodes;
        for (std::size_t i=1; i<=n_nodes; i++)
            nodes.push_back(new MeshLib::Node(static_cast<double>(i), static_cast<double>(i), i*0.5));

        std::vector<MeshLib::Element*> elements;
        for (std::size_t i=0; i<n_nodes-1; i++)
            elements.push_back(new MeshLib::Line(std::array<MeshLib::Node*,2>{{nodes[i], nodes[i+1]}}));

        _mesh = std::make_unique<MeshLib::Mesh>("TestMesh", nodes, elements);
    }

protected:
    std::unique_ptr<MeshLib::Mesh> _mesh;
};
// Project to parallels of XY plane
TEST_F(ProjectionTest, ProjectToXY)
{
    MathLib::Vector3 normal (0,0,1);
    std::size_t const n_nodes (_mesh->getNumberOfNodes());
    for (std::size_t p=0; p<10; p++)
    {
        MathLib::Point3d origin (std::array<double,3>{{0,0,static_cast<double>(p)}});
        MeshLib::Mesh* result = MeshLib::projectMeshOntoPlane(*_mesh, origin, normal);
        for (std::size_t i=0; i<n_nodes; i++)
            ASSERT_NEAR(static_cast<double>(p), (*result->getNode(i))[2], std::numeric_limits<double>::epsilon());
        delete result;
    }
}

// Project to parallels of XZ plane
TEST_F(ProjectionTest, ProjectToXZ)
{
    MathLib::Vector3 normal (0,1,0);
    std::size_t const n_nodes (_mesh->getNumberOfNodes());
    for (std::size_t p=0; p<10; p++)
    {
        MathLib::Point3d origin (std::array<double,3>{{0,static_cast<double>(p),0}});
        MeshLib::Mesh* result = MeshLib::projectMeshOntoPlane(*_mesh, origin, normal);
        for (std::size_t i=0; i<n_nodes; i++)
            ASSERT_NEAR(static_cast<double>(p), (*result->getNode(i))[1], std::numeric_limits<double>::epsilon());
        delete result;
    }
}

// Project to parallels of YZ plane
TEST_F(ProjectionTest, ProjectToYZ)
{
    MathLib::Vector3 normal (1,0,0);
    std::size_t const n_nodes (_mesh->getNumberOfNodes());
    for (std::size_t p=0; p<10; p++)
    {
        MathLib::Point3d origin (std::array<double,3>{{static_cast<double>(p),0,0}});
        MeshLib::Mesh* result = MeshLib::projectMeshOntoPlane(*_mesh, origin, normal);
        for (std::size_t i=0; i<n_nodes; i++)
            ASSERT_NEAR(static_cast<double>(p), (*result->getNode(i))[0], std::numeric_limits<double>::epsilon());
        delete result;
    }
}

// Sign of normal vector does not matter.
TEST_F(ProjectionTest, NormalDirection)
{
    MathLib::Vector3 normal_p (0,0,1);
    MathLib::Vector3 normal_n (0,0,-1);
    std::size_t const n_nodes (_mesh->getNumberOfNodes());
    MathLib::Point3d origin (std::array<double,3>{{0,0,0}});
    MeshLib::Mesh* result_p = MeshLib::projectMeshOntoPlane(*_mesh, origin, normal_p);
    MeshLib::Mesh* result_n = MeshLib::projectMeshOntoPlane(*_mesh, origin, normal_n);
    for (std::size_t i=0; i<n_nodes; i++)
        ASSERT_EQ((*result_p->getNode(i))[2], (*result_n->getNode(i))[2]);
    delete result_p;
    delete result_n;
}

// Length of normal does not matter (it's normalised within the method)
TEST_F(ProjectionTest, NormalLength)
{
    MathLib::Point3d origin (std::array<double,3>{{0,0,0}});
    MathLib::Vector3 normal (0,0,1);
    std::size_t const n_nodes (_mesh->getNumberOfNodes());
    MeshLib::Mesh* result = MeshLib::projectMeshOntoPlane(*_mesh, origin, normal);
    for (std::size_t p=2; p<10; p++)
    {
        normal[2] = static_cast<double>(p);
        MeshLib::Mesh* result_p = MeshLib::projectMeshOntoPlane(*_mesh, origin, normal);
        for (std::size_t i=0; i<n_nodes; i++)
            ASSERT_EQ((*result->getNode(i))[2], (*result_p->getNode(i))[2]);
        delete result_p;
    }
    delete result;
}
