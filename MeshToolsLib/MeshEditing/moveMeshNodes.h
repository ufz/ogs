/**
 * \file
 * \author Thomas Fischer
 * \date Feb 03, 2014
 * \brief Functionality to move mesh nodes using a given displacement vec.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
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
