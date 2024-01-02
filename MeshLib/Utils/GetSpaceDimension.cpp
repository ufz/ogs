/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 8, 2021, 12:16 PM
 */

#include "GetSpaceDimension.h"

#include <Eigen/Core>
#include <algorithm>
#include <limits>

#include "BaseLib/Error.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
int getSpaceDimension(std::vector<Node*> const& nodes)
{
    Eigen::Vector3d x_magnitude = Eigen::Vector3d::Zero();
    auto const x_ref = nodes[0]->asEigenVector3d();

    for (auto const& node : nodes)
    {
        auto const x = node->asEigenVector3d();
        x_magnitude += (x - x_ref).cwiseAbs();
    }

    // Z coordinate norm is not zero whichever 1D, 2D or 3D mesh:
    if (x_magnitude[2] > 0)
    {
        return 3;
    }

    const auto computed_dimension = std::count_if(
        x_magnitude.begin(), x_magnitude.end(),
        [](const double x_i_magnitude)
        { return x_i_magnitude > std::numeric_limits<double>::epsilon(); });

    // 1D mesh in y direction:
    if (computed_dimension == 1 && x_magnitude[1] > 0)
    {
        return 2;
    }

    return static_cast<int>(computed_dimension);
}
};  // namespace MeshLib
