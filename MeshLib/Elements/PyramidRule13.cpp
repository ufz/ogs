/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PyramidRule13.h"

#include "BaseLib/Logging.h"
#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib
{
const unsigned PyramidRule13::n_face_nodes[5] = {6, 6, 6, 6, 8};

const Element* PyramidRule13::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        unsigned nFaceNodes(PyramidRule13::n_face_nodes[i]);
        auto** nodes = new Node*[nFaceNodes];
        for (unsigned j = 0; j < nFaceNodes; j++)
        {
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        }

        if (i < 4)
        {
            return new Tri6(nodes, e->getID());
        }

        return new Quad8(nodes, e->getID());
    }
    ERR("Error in MeshLib::Element::getFace() - Index {:d} does not exist.", i);
    return nullptr;
}

}  // end namespace MeshLib
