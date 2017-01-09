/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "QuadRule8.h"

namespace MeshLib {

const unsigned QuadRule8::edge_nodes[4][3] =
{
        {0, 1, 4}, // Edge 0
        {1, 2, 5}, // Edge 1
        {2, 3, 6}, // Edge 2
        {0, 3, 7}  // Edge 3
};

} // end namespace MeshLib
