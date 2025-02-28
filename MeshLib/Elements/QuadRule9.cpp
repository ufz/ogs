/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "QuadRule9.h"

namespace MeshLib
{
const unsigned QuadRule9::edge_nodes[4][3] = {
    {0, 1, 4},  // Edge 0
    {1, 2, 5},  // Edge 1
    {2, 3, 6},  // Edge 2
    {3, 0, 7}   // Edge 3
};

}  // end namespace MeshLib
