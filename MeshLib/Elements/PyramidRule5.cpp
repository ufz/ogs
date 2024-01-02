/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PyramidRule5.h"

#include "BaseLib/Logging.h"
#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib
{
const unsigned PyramidRule5::face_nodes[5][4] = {
    {0, 1, 4, 99},  // Face 0
    {1, 2, 4, 99},  // Face 1
    {2, 3, 4, 99},  // Face 2
    {3, 0, 4, 99},  // Face 3
    {0, 3, 2, 1}    // Face 4
};

const unsigned PyramidRule5::edge_nodes[8][2] = {
    {0, 1},  // Edge 0
    {1, 2},  // Edge 1
    {2, 3},  // Edge 2
    {0, 3},  // Edge 3
    {0, 4},  // Edge 4
    {1, 4},  // Edge 5
    {2, 4},  // Edge 6
    {3, 4}   // Edge 7
};

const unsigned PyramidRule5::n_face_nodes[5] = {3, 3, 3, 3, 4};

const Element* PyramidRule5::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        unsigned nFaceNodes(PyramidRule5::n_face_nodes[i]);
        auto** nodes = new Node*[nFaceNodes];
        for (unsigned j = 0; j < nFaceNodes; j++)
        {
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        }

        if (i < 4)
        {
            return new Tri(nodes, e->getID());
        }

        return new Quad(nodes, e->getID());
    }
    ERR("Error in MeshLib::Element::getFace() - Index {:d} does not exist.", i);
    return nullptr;
}

}  // end namespace MeshLib
