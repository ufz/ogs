/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HexRule8.h"

#include <array>

#include "BaseLib/Logging.h"
#include "Line.h"
#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"
#include "Quad.h"

namespace MeshLib
{
const unsigned HexRule8::face_nodes[6][4] = {
    {0, 3, 2, 1},  // Face 0
    {0, 1, 5, 4},  // Face 1
    {1, 2, 6, 5},  // Face 2
    {2, 3, 7, 6},  // Face 3
    {3, 0, 4, 7},  // Face 4
    {4, 5, 6, 7}   // Face 5
};

const unsigned HexRule8::edge_nodes[12][2] = {
    {0, 1},  // Edge 0
    {1, 2},  // Edge 1
    {2, 3},  // Edge 2
    {0, 3},  // Edge 3
    {4, 5},  // Edge 4
    {5, 6},  // Edge 5
    {6, 7},  // Edge 6
    {4, 7},  // Edge 7
    {0, 4},  // Edge 8
    {1, 5},  // Edge 9
    {2, 6},  // Edge 10
    {3, 7}   // Edge 11
};

const Element* HexRule8::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        std::array<Node*, 4> nodes{};
        for (unsigned j = 0; j < 4; j++)
        {
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        }
        return new Quad(nodes, e->getID());
    }
    ERR("Error in MeshLib::Element::getFace() - Index {:d} does not exist.", i);
    return nullptr;
}
}  // end namespace MeshLib
