/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TriRule3.h"

namespace MeshLib
{
const unsigned TriRule3::edge_nodes[3][2] = {
    {0, 1},  // Edge 0
    {1, 2},  // Edge 1
    {2, 0},  // Edge 2
};
}  // end namespace MeshLib
