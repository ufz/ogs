/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class ISpatialFunction;

/**
 * Generate distributed node values from a given function
 *
 * @param func            a spatial function object
 * @param msh             a mesh object
 * @param vec_node_ids    a vector of mesh node ids where the function is evaluated
 * @return a vector of nodal values
 */
std::vector<double> generateNodeValueDistribution(
    const NumLib::ISpatialFunction &func,
    const MeshLib::Mesh &msh,
    const std::vector<std::size_t> &vec_node_ids);

} // NumLib
