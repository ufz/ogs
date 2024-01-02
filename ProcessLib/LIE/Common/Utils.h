/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{
/// compute physical coordinates from the given shape vector, i.e. from the
/// natural coordinates
template <typename Derived>
Eigen::Vector3d computePhysicalCoordinates(
    MeshLib::Element const& e, Eigen::MatrixBase<Derived> const& shape)
{
    Eigen::Vector3d pt = Eigen::Vector3d::Zero();
    for (unsigned i = 0; i < e.getNumberOfNodes(); i++)
    {
        MeshLib::Node const& node = *e.getNode(i);
        for (unsigned j = 0; j < 3; j++)
        {
            pt[j] += shape[i] * node[j];
        }
    }
    return pt;
}

}  // namespace LIE
}  // namespace ProcessLib
