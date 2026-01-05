// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "HexRule8.h"

#include <array>

#include "BaseLib/Logging.h"
#include "Line.h"
#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"
#include "Quad.h"

namespace MeshLib
{
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
