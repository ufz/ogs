/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 18, 2021, 12:31 PM
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
#include "MeshLib/Utils/GetElementRotationMatrices.h"
#include "MeshLib/Utils/GetSpaceDimension.h"

TEST(MeshLib, GetElementRotationMatrices3DMesh)
{
    // The memory of the nodes are allocated by new operator, and it is released
    // in the destructor of MeshLib::Mesh.

    // Construct a 3D mesh, which contains two inclined 2D triangle elements
    // and an inclined line element.
    std::vector<MeshLib::Node*> nodes(10);
    nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
    nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0, 1);
    nodes[2] = new MeshLib::Node(1.0, 1.0, 0.0, 2);
    nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0, 3);

    nodes[4] = new MeshLib::Node(0.0, 0.0, 1.0, 4);
    nodes[5] = new MeshLib::Node(1.0, 0.0, 1.0, 5);
    nodes[6] = new MeshLib::Node(1.0, 1.0, 1.0, 6);
    nodes[7] = new MeshLib::Node(0.0, 1.0, 1.0, 7);

    nodes[8] = new MeshLib::Node(0.0, 0.5, 1.5, 8);
    nodes[9] = new MeshLib::Node(1.0, 0.5, 1.5, 9);

    std::vector<MeshLib::Element*> elements;

    // One hexahedral element:
    std::array<MeshLib::Node*, 8> hex_element_nodes;
    std::copy_n(nodes.begin(), 8, hex_element_nodes.begin());
    elements.push_back(new MeshLib::Hex(hex_element_nodes));

    // Two inclined triangle elements:
    std::array<MeshLib::Node*, 3> tri_element_nodes{nodes[6], nodes[7],
                                                    nodes[8]};
    elements.push_back(new MeshLib::Tri(tri_element_nodes));
    tri_element_nodes[0] = nodes[6];
    tri_element_nodes[1] = nodes[8];
    tri_element_nodes[2] = nodes[9];
    elements.push_back(new MeshLib::Tri(tri_element_nodes));

    // One inclined line element:
    std::array<MeshLib::Node*, 2> line_element_nodes{nodes[6], nodes[8]};
    elements.push_back(new MeshLib::Line(line_element_nodes));

    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    meshes.push_back(
        std::make_unique<MeshLib::Mesh>("a_mesh", nodes, elements));

    int const space_dimension = MeshLib::getSpaceDimension(nodes);
    auto const element_rotation_matrices = MeshLib::getElementRotationMatrices(
        space_dimension, meshes[0]->getDimension(), meshes[0]->getElements());

    // First element, the hexahedral element has identity matrix:
    EXPECT_EQ(9, element_rotation_matrices[0].size());

    // Second and third elements, the inclined triangle elements on the
    // same plane:
    auto const& rotation_matrix_tri1 = element_rotation_matrices[1];
    auto const& rotation_matrix_tri2 = element_rotation_matrices[2];

    double const diff = (rotation_matrix_tri1 - rotation_matrix_tri2).norm();
    ASSERT_LE(diff, 1e-10);

    Eigen::VectorXd b(3);
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 1.0;

    // Projection test
    Eigen::VectorXd const b_local = rotation_matrix_tri1.transpose() * b;
    EXPECT_EQ(elements[1]->getDimension(), b_local.size());
    EXPECT_EQ(0.0, b_local[0]);
    double const expected_b_local1 = -std::sqrt(2.0) / 2.0;
    EXPECT_LE(std::fabs(b_local[1] - expected_b_local1), 1.e-16);

    // Forth element, the inclined line element:
    Eigen::VectorXd const b_local_1D =
        element_rotation_matrices[3].transpose() * b;
    EXPECT_EQ(elements[3]->getDimension(), b_local_1D.size());

    double const* const x_6 = nodes[6]->getCoords();
    double const* const x_8 = nodes[8]->getCoords();
    // b_local_1D = |b| (x_8-x_6) * b/(|(x_8-x_6)| |b|)
    double dx[3];
    for (int i = 0; i < 3; i++)
    {
        dx[i] = x_8[i] - x_6[i];
    }

    double const dx_norm =
        std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    double const expected_b_local_1D = dx[2] / dx_norm;
    EXPECT_LE(std::fabs(expected_b_local_1D - b_local_1D[0]), 1.e-16);

    // Test of rotation of a local vector to the global system
    Eigen::VectorXd const local_1D_to_global =
        element_rotation_matrices[3] * b_local_1D;
    EXPECT_EQ(meshes[0]->getDimension(), local_1D_to_global.size());

    double const local_1D_to_global_norm =
        std::sqrt(local_1D_to_global[0] * local_1D_to_global[0] +
                  local_1D_to_global[1] * local_1D_to_global[1] +
                  local_1D_to_global[2] * local_1D_to_global[2]);
    for (int i = 0; i < 3; i++)
    {
        EXPECT_LE(std::fabs(dx[i] / dx_norm -
                            local_1D_to_global[i] / local_1D_to_global_norm),
                  1.e-15);
    }
}
