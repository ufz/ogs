/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TetRule4.h"

#include <array>

#include "BaseLib/Logging.h"
#include "MeshLib/Node.h"
#include "Tri.h"

namespace MeshLib
{
const unsigned TetRule4::face_nodes[4][3] = {
    {0, 2, 1},  // Face 0
    {0, 1, 3},  // Face 1
    {1, 2, 3},  // Face 2
    {2, 0, 3}   // Face 3
};

const unsigned TetRule4::edge_nodes[6][2] = {
    {0, 1},  // Edge 0
    {1, 2},  // Edge 1
    {0, 2},  // Edge 2
    {0, 3},  // Edge 3
    {1, 3},  // Edge 4
    {2, 3}   // Edge 5
};

const Element* TetRule4::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        std::array<Node*, 3> nodes{};
        for (unsigned j = 0; j < 3; j++)
        {
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        }
        return new Tri(nodes, e->getID());
    }
    ERR("Error in MeshLib::Element::getFace() - Index {:d} does not exist.", i);
    return nullptr;
}
}  // end namespace MeshLib
