// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PyramidRule5.h"

#include "BaseLib/Logging.h"
#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib
{
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
