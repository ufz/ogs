/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Distribution.h"

#include <algorithm>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "NumLib/Function/ISpatialFunction.h"

namespace NumLib
{

std::vector<double> generateNodeValueDistribution(
    const NumLib::ISpatialFunction &func,
    const MeshLib::Mesh &msh,
    const std::vector<std::size_t> &vec_node_ids)
{
    // evaluate a given function with nodal coordinates
    std::vector<double> vec_values(vec_node_ids.size());
    std::transform(vec_node_ids.begin(), vec_node_ids.end(), vec_values.begin(),
        [&func, &msh](std::size_t node_id)
        {
            return func(*msh.getNode(node_id));
        });
    return vec_values;
}

} //NumLib
