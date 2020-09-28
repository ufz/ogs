/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 5, 2021, 12:46 PM
 */

#include "SetMeshSpaceDimension.h"

#include <algorithm>
#include <array>
#include <limits>

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
unsigned getSpaceDimension(Mesh const& mesh)
{
    std::array x_magnitude = {0.0, 0.0, 0.0};

    auto const nodes = mesh.getNodes();
    for (auto const& node : nodes)
    {
        auto const x = node->getCoords();
        for (int i = 0; i < 3; i++)
        {
            x_magnitude[i] += std::fabs(x[i]);
        }
    }

    return static_cast<unsigned>(std::count_if(
        x_magnitude.begin(), x_magnitude.end(), [](const double x_i_magnitude) {
            return x_i_magnitude > std::numeric_limits<double>::epsilon();
        }));
}

void setMeshSpaceDimension(std::vector<std::unique_ptr<Mesh>> const& meshes)
{
    // Get the space dimension from the bulk mesh:
    auto const space_dimension = getSpaceDimension(*meshes[0]);
    for (auto& mesh : meshes)
    {
        auto elements = mesh->getElements();
        for (auto element : elements)
        {
            element->space_dimension_ = space_dimension;
        }
    }
}
};  // namespace MeshLib
