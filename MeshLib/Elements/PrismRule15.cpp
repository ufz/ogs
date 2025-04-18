/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PrismRule15.h"

#include "BaseLib/Logging.h"
#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib
{
const unsigned PrismRule15::face_nodes[5][8] = {
    {0, 2, 1, 8, 7, 6, 99, 99},   // Face 0
    {0, 1, 4, 3, 6, 13, 9, 12},   // Face 1
    {1, 2, 5, 4, 7, 14, 10, 13},  // Face 2
    {2, 0, 3, 5, 8, 12, 11, 14},  // Face 3
    {3, 4, 5, 9, 10, 11, 99, 99}  // Face 4
};

const unsigned PrismRule15::edge_nodes[9][3] = {
    {0, 1, 6},   // Edge 0
    {1, 2, 7},   // Edge 1
    {0, 2, 8},   // Edge 2
    {0, 3, 12},  // Edge 3
    {1, 4, 13},  // Edge 4
    {2, 5, 14},  // Edge 5
    {3, 4, 9},   // Edge 6
    {4, 5, 10},  // Edge 7
    {3, 5, 11}   // Edge 8
};

const unsigned PrismRule15::n_face_nodes[5] = {6, 8, 8, 8, 6};

const Element* PrismRule15::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        unsigned nFaceNodes(PrismRule15::n_face_nodes[i]);
        auto** nodes = new Node*[nFaceNodes];
        for (unsigned j = 0; j < nFaceNodes; j++)
        {
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        }

        if (i == 0 || i == 4)
        {
            return new Tri6(nodes, e->getID());
        }

        return new Quad8(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index {:d} does not exist.", i);
    return nullptr;
}

}  // end namespace MeshLib
