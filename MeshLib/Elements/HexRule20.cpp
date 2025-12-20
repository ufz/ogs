// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "HexRule20.h"

#include <array>

#include "BaseLib/Logging.h"
#include "Line.h"
#include "MeshLib/Node.h"
#include "Quad.h"

namespace MeshLib
{

const Element* HexRule20::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        std::array<Node*, 8> nodes{};
        for (unsigned j = 0; j < 8; j++)
        {
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        }
        return new Quad8(nodes, e->getID());
    }
    ERR("Error in MeshLib::Element::getFace() - Index {:d} does not exist.", i);
    return nullptr;
}

}  // end namespace MeshLib
