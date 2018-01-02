/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PyramidRule13.h"

#include <logog/include/logog.hpp>

#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib {

const unsigned PyramidRule13::face_nodes[5][8] =
{
    {0, 1, 4, 5, 10, 9, 99, 99}, // Face 0
    {1, 2, 4, 6, 11, 10, 99, 99}, // Face 1
    {2, 3, 4, 7, 12, 11, 99, 99}, // Face 2
    {3, 0, 4, 8, 9, 12, 99, 99}, // Face 3
    {0, 3, 2, 1, 8, 7, 6, 5}  // Face 4
};

const unsigned PyramidRule13::edge_nodes[8][3] =
{
    {0, 1, 5}, // Edge 0
    {1, 2, 6}, // Edge 1
    {2, 3, 7}, // Edge 2
    {0, 3, 8}, // Edge 3
    {0, 4, 9}, // Edge 4
    {1, 4, 10}, // Edge 5
    {2, 4, 11}, // Edge 6
    {3, 4, 12}  // Edge 7
};

const unsigned PyramidRule13::n_face_nodes[5] = { 6, 6, 6, 6, 8 };

const Element* PyramidRule13::getFace(const Element* e, unsigned i)
{
    if (i<e->getNumberOfFaces())
    {
        unsigned nFaceNodes(PyramidRule13::n_face_nodes[i]);
        auto** nodes = new Node*[nFaceNodes];
        for (unsigned j=0; j<nFaceNodes; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));

        if (i < 4)
            return new Tri6(nodes);

        return new Quad8(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

} // end namespace MeshLib
