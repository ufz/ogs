/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TetRule10.h"

#include <array>

#include <logog/include/logog.hpp>

#include "MeshLib/Node.h"
#include "Tri.h"

namespace MeshLib
{

const unsigned TetRule10::face_nodes[4][6] =
{
    {0, 2, 1, 6, 5, 4}, // Face 0
    {0, 1, 3, 4, 8, 7}, // Face 1
    {1, 2, 3, 5, 9, 8}, // Face 2
    {2, 0, 3, 6, 7, 9}  // Face 3
};

const unsigned TetRule10::edge_nodes[6][3] =
{
    {0, 1, 4}, // Edge 0
    {1, 2, 5}, // Edge 1
    {0, 2, 6}, // Edge 2
    {0, 3, 7}, // Edge 3
    {1, 3, 8}, // Edge 4
    {2, 3, 9}  // Edge 5
};

const Element* TetRule10::getFace(const Element* e, unsigned i)
{
    if (i<n_faces)
    {
        std::array<Node*,6> nodes;
        for (unsigned j=0; j<6; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        return new Tri6(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

} // end namespace MeshLib
