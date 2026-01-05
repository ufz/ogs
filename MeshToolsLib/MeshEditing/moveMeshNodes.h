// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

#include "MathLib/Point3d.h"

namespace MeshToolsLib
{

/**
 * Function that moves mesh nodes.
 *
 * The function iterates over all mesh nodes between the
 * begin and the end iterator and moves them using the
 * given displacement.
 * \param begin begin iterator
 * \param end end iterator
 * \param displacement the displacement to use
 */
template <typename Iterator>
void moveMeshNodes(Iterator begin,
                   Iterator end,
                   Eigen::Vector3d const& displacement)
{
    std::for_each(begin, end,
                  [&displacement](MathLib::Point3d* node)
                  { node->asEigenVector3d() += displacement; });
};

}  // namespace MeshToolsLib
