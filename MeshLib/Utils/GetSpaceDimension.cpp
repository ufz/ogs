/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 8, 2021, 12:16 PM
 */

#include "GetSpaceDimension.h"

#include <algorithm>
#include <array>
#include <limits>

#include "BaseLib/Error.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
int getSpaceDimension(std::vector<Node*> const& nodes)
{
    std::array x_magnitude = {0.0, 0.0, 0.0};

    double const* const x_ref = nodes[0]->getCoords();
    for (auto const& node : nodes)
    {
        auto const x = node->getCoords();
        for (int i = 0; i < 3; i++)
        {
            x_magnitude[i] += std::fabs(x[i] - x_ref[i]);
        }
    }

    return static_cast<int>(std::count_if(
        x_magnitude.begin(), x_magnitude.end(),
        [](const double x_i_magnitude)
        { return x_i_magnitude > std::numeric_limits<double>::epsilon(); }));
}
};  // namespace MeshLib
