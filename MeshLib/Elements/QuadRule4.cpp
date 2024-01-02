/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "QuadRule4.h"

#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
const unsigned QuadRule4::edge_nodes[4][2] = {
    {0, 1},  // Edge 0
    {1, 2},  // Edge 1
    {2, 3},  // Edge 2
    {3, 0}   // Edge 3
};
}  // end namespace MeshLib
