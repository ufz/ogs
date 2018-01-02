/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HexRule20.h"

#include <array>

#include <logog/include/logog.hpp>

#include "MeshLib/Node.h"
#include "Quad.h"
#include "Line.h"

namespace MeshLib {

const unsigned HexRule20::face_nodes[6][8] =
{
    {0, 3, 2, 1, 11, 10,  9,  8}, // Face 0
    {0, 1, 5, 4,  8, 17, 12, 11}, // Face 1
    {1, 2, 6, 5,  9, 18, 13, 17}, // Face 2
    {2, 3, 7, 6, 10, 19, 14, 18}, // Face 3
    {3, 0, 4, 7, 11, 16, 15, 19}, // Face 4
    {4, 5, 6, 7, 12, 13, 14, 15}  // Face 5
};

const unsigned HexRule20::edge_nodes[12][3] =
{
    {0, 1, 8}, // Edge 0
    {1, 2, 9}, // Edge 1
    {2, 3, 10}, // Edge 2
    {0, 3, 11}, // Edge 3
    {4, 5, 12}, // Edge 4
    {5, 6, 13}, // Edge 5
    {6, 7, 14}, // Edge 6
    {4, 7, 15}, // Edge 7
    {0, 4, 16}, // Edge 8
    {1, 5, 17}, // Edge 9
    {2, 6, 18}, // Edge 10
    {3, 7, 19}  // Edge 11
};

const Element* HexRule20::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        std::array<Node*, 8> nodes;
        for (unsigned j=0; j<8; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        return new Quad8(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

} // end namespace MeshLib
