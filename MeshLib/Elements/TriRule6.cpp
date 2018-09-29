/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TriRule6.h"

namespace MeshLib {

const unsigned TriRule6::edge_nodes[3][3] =
{
        {0, 1, 3}, // Edge 0
        {1, 2, 4}, // Edge 1
        {2, 0, 5}, // Edge 2
};

} // end namespace MeshLib
