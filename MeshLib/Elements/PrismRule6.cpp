/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PrismRule6.h"

#include "BaseLib/Logging.h"
#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib
{
const unsigned PrismRule6::face_nodes[5][4] = {
    {0, 2, 1, 99},  // Face 0
    {0, 1, 4, 3},   // Face 1
    {1, 2, 5, 4},   // Face 2
    {2, 0, 3, 5},   // Face 3
    {3, 4, 5, 99}   // Face 4
};

const unsigned PrismRule6::edge_nodes[9][2] = {
    {0, 1},  // Edge 0
    {1, 2},  // Edge 1
    {0, 2},  // Edge 2
    {0, 3},  // Edge 3
    {1, 4},  // Edge 4
    {2, 5},  // Edge 5
    {3, 4},  // Edge 6
    {4, 5},  // Edge 7
    {3, 5}   // Edge 8
};

const unsigned PrismRule6::n_face_nodes[5] = {3, 4, 4, 4, 3};

const Element* PrismRule6::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        unsigned nFaceNodes(PrismRule6::n_face_nodes[i]);
        auto** nodes = new Node*[nFaceNodes];
        for (unsigned j = 0; j < nFaceNodes; j++)
        {
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        }

        if (i == 0 || i == 4)
        {
            return new Tri(nodes, e->getID());
        }

        return new Quad(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index {:d} does not exist.", i);
    return nullptr;
}
}  // end namespace MeshLib
