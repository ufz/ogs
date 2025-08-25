/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TetRule10.h"

#include <array>

#include "BaseLib/Logging.h"
#include "MeshLib/Node.h"
#include "Tri.h"

namespace MeshLib
{
const Element* TetRule10::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        std::array<Node*, 6> nodes{};
        for (unsigned j = 0; j < 6; j++)
        {
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        }
        return new Tri6(nodes, e->getID());
    }
    ERR("Error in MeshLib::Element::getFace() - Index {:d} does not exist.", i);
    return nullptr;
}

}  // end namespace MeshLib
