/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <vector>

#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/moveMeshNodes.h"

TEST(MeshLib, moveMeshNodes)
{
    /* initialize random seed: */
    srand(static_cast<unsigned>(time(nullptr)));

    std::size_t const size(16384);

    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Node*> nodes_copy;
    nodes.resize(size);
    nodes_copy.resize(size);

    /* put nodes with random coords into vectors */
    for (std::size_t k(0); k < size; k++)
    {
        nodes[k] = new MeshLib::Node(rand(), rand(), rand());
        nodes_copy[k] = new MeshLib::Node(*nodes[k]);
    }

    /* create random displacement */
    Eigen::Vector3d displacement{static_cast<double>(rand()),
                                 static_cast<double>(rand()),
                                 static_cast<double>(rand())};

    /* move the mesh node */
    MeshToolsLib::moveMeshNodes(nodes.begin(), nodes.end(), displacement);

    /* reverse the direction of displacement */
    displacement *= -1.0;

    /* move the mesh node back */
    MeshToolsLib::moveMeshNodes(nodes.begin(), nodes.end(), displacement);

    /* check the result */
    double const eps(std::numeric_limits<double>::epsilon());
    for (std::size_t k(0); k < size; k++)
    {
        EXPECT_NEAR((*nodes_copy[0])[0], (*nodes[0])[0], eps);
        EXPECT_NEAR((*nodes_copy[0])[1], (*nodes[0])[1], eps);
        EXPECT_NEAR((*nodes_copy[0])[2], (*nodes[0])[2], eps);
    }

    for (auto n : nodes)
    {
        delete n;
    }
    for (auto n : nodes_copy)
    {
        delete n;
    }
}
