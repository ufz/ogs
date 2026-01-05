// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "TetRule4.h"

#include <array>

#include "BaseLib/Logging.h"
#include "MeshLib/Node.h"
#include "Tri.h"

namespace MeshLib
{
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
