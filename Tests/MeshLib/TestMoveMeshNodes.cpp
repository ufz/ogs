/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <cstdlib>
#include <ctime>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "MeshLib/MeshEditing/moveMeshNodes.h"
#include "MeshLib/Node.h"

TEST(MeshLib, moveMeshNodes)
{
    /* initialize random seed: */
    srand(static_cast<unsigned>(time(nullptr)));

    std::size_t const size (16384);

    std::vector<MeshLib::Node*> nodes, nodes_copy;
    nodes.resize(size);
    nodes_copy.resize(size);

    /* put nodes with random coords into vectors */
    for (std::size_t k(0); k<size; k++) {
        nodes[k] = new MeshLib::Node(rand(), rand(), rand());
        nodes_copy[k] = new MeshLib::Node(* nodes[k]);
    }

    /* create random displacement */
    MeshLib::Node displacement(rand(), rand(), rand());

    /* move the mesh node */
    MeshLib::moveMeshNodes(nodes.begin(), nodes.end(), displacement);

    /* reverse the direction of displacement */
    displacement[0] *= -1.0;
    displacement[1] *= -1.0;
    displacement[2] *= -1.0;

    /* move the mesh node back */
    MeshLib::moveMeshNodes(nodes.begin(), nodes.end(), displacement);

    /* check the result */
    double const eps(std::numeric_limits<double>::epsilon());
    for (std::size_t k(0); k<size; k++) {
        EXPECT_NEAR((*nodes_copy[0])[0], (*nodes[0])[0], eps);
        EXPECT_NEAR((*nodes_copy[0])[1], (*nodes[0])[1], eps);
        EXPECT_NEAR((*nodes_copy[0])[2], (*nodes[0])[2], eps);
    }

    for (auto n : nodes)
        delete n;
    for (auto n : nodes_copy)
        delete n;
}

