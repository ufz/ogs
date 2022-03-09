/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 10, 2021, 3:09 PM
 */

#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/Is2DMeshOnRotatedVerticalPlane.h"
#include "MeshLib/Utils/SetMeshSpaceDimension.h"
template <typename ElementType, int NodeNumber>
void testAsNoVerticalPlaneMesh(
    std::string const& error_message,
    const unsigned expected_space_dimension,
    std::array<MeshLib::Node*, NodeNumber> const& element_nodes)
{
    std::vector<MeshLib::Element*> elements{new ElementType(element_nodes)};

    std::vector<MeshLib::Node*> nodes(element_nodes.begin(),
                                      element_nodes.end());

    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    meshes.push_back(
        std::make_unique<MeshLib::Mesh>("a_mesh", nodes, elements));

    MeshLib::setMeshSpaceDimension(meshes);
    ASSERT_EQ(expected_space_dimension, elements[0]->space_dimension_);

    try
    {
        MeshLib::is2DMeshOnRotatedVerticalPlane(*meshes[0]);
    }
    catch (std::exception& e)
    {
        EXPECT_EQ(e.what(), error_message);
    }
}

void testAsVerticalPlaneMesh(const unsigned expected_space_dimension,
                             std::array<MeshLib::Node*, 4> const& element_nodes)
{
    std::vector<MeshLib::Element*> elements{new MeshLib::Quad(element_nodes)};

    std::vector<MeshLib::Node*> nodes(element_nodes.begin(),
                                      element_nodes.end());

    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    meshes.push_back(
        std::make_unique<MeshLib::Mesh>("a_mesh", nodes, elements));

    MeshLib::setMeshSpaceDimension(meshes);

    if (expected_space_dimension == 3u)
    {
        ASSERT_EQ(3u, elements[0]->space_dimension_);

        ASSERT_TRUE(MeshLib::is2DMeshOnRotatedVerticalPlane(*meshes[0]));
    }
    else
    {
        ASSERT_EQ(2u, elements[0]->space_dimension_);

        ASSERT_FALSE(MeshLib::is2DMeshOnRotatedVerticalPlane(*meshes[0]));
    }
}

TEST(MeshLib, Is2DMeshOnRotatedVerticalPlane)
{
    // The memory of the nodes are allocated by new operator, and it is released
    // in the destructor of MeshLib::Mesh.

    // 3D mesh:
    {
        std::array<MeshLib::Node*, 4> element_nodes;
        element_nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
        element_nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0, 1);
        element_nodes[2] = new MeshLib::Node(0.0, 1.0, 0.0, 2);
        element_nodes[3] = new MeshLib::Node(0.0, 0.0, 1.0, 3);

        std::string const error_message =
            "A 2D mesh is required for this "
            "computation but the provided mesh, "
            "mesh a_mesh, has 3D elements.";

        auto const expected_space_dimension = 3u;
        testAsNoVerticalPlaneMesh<MeshLib::Tet, 4>(
            error_message, expected_space_dimension, element_nodes);
    }

    // Inclined 2D mesh but it is not on a plane rotated around the vertical
    // axis:
    {
        std::array<MeshLib::Node*, 4> element_nodes;
        element_nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
        element_nodes[1] = new MeshLib::Node(1.0, 1.0, 0.0, 1);
        element_nodes[2] = new MeshLib::Node(0.2, 0.5, 1.0, 2);
        element_nodes[3] = new MeshLib::Node(1.2, 1.5, 0.0, 3);

        std::string const error_message =
            "2D Mesh a_mesh is on an inclined plane, which is neither a "
            "vertical nor horizontal plane that is required for the present "
            "computation.";

        auto const expected_space_dimension = 3u;
        testAsNoVerticalPlaneMesh<MeshLib::Quad, 4>(
            error_message, expected_space_dimension, element_nodes);
    }

    // 1D mesh:
    {
        std::array<MeshLib::Node*, 2> element_nodes;
        element_nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
        element_nodes[1] = new MeshLib::Node(0.2, 0.5, 1.0, 1);

        std::string const error_message =
            "A 2D mesh is required for this "
            "computation but the provided mesh, "
            "mesh a_mesh, has 1D elements.";

        auto const expected_space_dimension = 3u;
        testAsNoVerticalPlaneMesh<MeshLib::Line, 2>(
            error_message, expected_space_dimension, element_nodes);
    }

    // 2D mesh on x-y plane:
    {
        std::array<MeshLib::Node*, 4> element_nodes;
        element_nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
        element_nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0, 1);
        element_nodes[2] = new MeshLib::Node(1.0, 1.0, 0.0, 2);
        element_nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0, 3);

        auto const expected_space_dimension = 2u;
        testAsVerticalPlaneMesh(expected_space_dimension, element_nodes);
    }

    // 2D mesh on x-z plane:
    {
        std::array<MeshLib::Node*, 4> element_nodes;
        element_nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
        element_nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0, 1);
        element_nodes[2] = new MeshLib::Node(1.0, 0.0, 1.0, 2);
        element_nodes[3] = new MeshLib::Node(0.0, 0.0, 1.0, 3);

        auto const expected_space_dimension = 3u;
        testAsVerticalPlaneMesh(expected_space_dimension, element_nodes);
    }

    // 2D mesh on the plane rotated around y axis:
    {
        std::array<MeshLib::Node*, 4> element_nodes;
        element_nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
        element_nodes[1] = new MeshLib::Node(0.2, 0.0, 0.5, 1);
        element_nodes[2] = new MeshLib::Node(0.2, 1.0, 0.5, 2);
        element_nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0, 3);

        auto const expected_space_dimension = 3u;
        testAsVerticalPlaneMesh(expected_space_dimension, element_nodes);
    }

    // 2D mesh on the plane rotated around z axis:
    {
        std::array<MeshLib::Node*, 4> element_nodes;
        element_nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
        element_nodes[1] = new MeshLib::Node(0.2, 0.5, 0.0, 1);
        element_nodes[2] = new MeshLib::Node(0.2, 0.5, 1.0, 2);
        element_nodes[3] = new MeshLib::Node(0.0, 0.0, 1.0, 3);

        auto const expected_space_dimension = 3u;
        testAsVerticalPlaneMesh(expected_space_dimension, element_nodes);
    }
}
