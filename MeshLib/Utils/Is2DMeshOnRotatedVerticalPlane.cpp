/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 22, 2021, 11:52 AM
 */

#include "Is2DMeshOnRotatedVerticalPlane.h"

#include <algorithm>
#include <limits>

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
bool is2DMeshOnRotatedVerticalPlane(Mesh const& mesh)
{
    auto const mesh_dimension = mesh.getDimension();
    if (mesh_dimension != 2)
    {
        OGS_FATAL(
            "A 2D mesh is required for this computation but the provided mesh, "
            "mesh {:s}, has {:d}D elements.",
            mesh.getName(), mesh_dimension);
    }

    auto const& elements = mesh.getElements();

    bool const has_inclined_element =
        std::any_of(elements.begin(), elements.end(),
                    [&mesh_dimension](auto const& element)
                    { return element->space_dimension_ != mesh_dimension; });

    if (!has_inclined_element)
    {
        return false;
    }

    const bool is_rotated_around_y_axis = std::all_of(
        elements.cbegin(), elements.cend(),
        [](auto const& element)
        {
            // 3 nodes are enough to make up a plane.
            auto const x1 = element->getNode(0)->data();
            auto const x2 = element->getNode(1)->data();
            auto const x3 = element->getNode(2)->data();

            double const a0 = x2[0] - x1[0];
            double const a2 = x2[2] - x1[2];

            double const b0 = x3[0] - x1[0];
            double const b2 = x3[2] - x1[2];

            double const e_n_1 = -a0 * b2 + a2 * b0;
            return std::fabs(e_n_1) < std::numeric_limits<double>::epsilon();
        });

    const bool is_rotated_around_z_axis = std::all_of(
        elements.cbegin(), elements.cend(),
        [](auto const& element)
        {
            // 3 nodes are enough to make up a plane.
            auto const x1 = element->getNode(0)->data();
            auto const x2 = element->getNode(1)->data();
            auto const x3 = element->getNode(2)->data();

            double const a0 = x2[0] - x1[0];
            double const a1 = x2[1] - x1[1];

            double const b0 = x3[0] - x1[0];
            double const b1 = x3[1] - x1[1];

            double const e_n_2 = a0 * b1 - a1 * b0;
            return std::fabs(e_n_2) < std::numeric_limits<double>::epsilon();
        });

    if (!(is_rotated_around_y_axis || is_rotated_around_z_axis))
    {
        OGS_FATAL(
            "2D Mesh {:s} is on an inclined plane, which is neither a vertical "
            "nor horizontal plane that is required for the present "
            "computation.",
            mesh.getName());
    }

    return true;
}
};  // namespace MeshLib
