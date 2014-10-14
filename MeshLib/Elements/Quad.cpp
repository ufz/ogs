/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Quad.h"

namespace MeshLib
{
namespace detail
{
constexpr unsigned QuadEdgeLinearNodes::_edge_nodes[4][2];
constexpr unsigned QuadEdgeQuadraticNodes::_edge_nodes[4][3];
}
}

