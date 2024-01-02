/*
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/// See Tensor type for details.
constexpr int tensorSize(int dim)
{
    if (dim == 1)
    {
        return 3;  // Diagonal entries.
    }
    if (dim == 2)
    {
        return 5;  // 2x2 matrix and the 3rd diagonal entry.
    }
    if (dim == 3)
    {
        return 9;  // Full 3x3 matrix.
    }
    OGS_FATAL("Tensor size for dimension {} is not defined.", dim);
}

/// The tensor's components in 3D case are ordered in the usual way:
/// (1,1), (1,2), (1,3)
/// (2,1), (2,2), (2,3)
/// (3,1), (3,2), (3,3).
///
/// For the 2D case the 2x2 block is as usual and is followed by the (3,3)
/// component:
/// (1,1), (1,2),
/// (2,1), (2,2),
///               (3,3).
///
/// For the 1D case only the diagonal is stored:
/// (1,1),
///        (2,2),
///               (3,3).
template <int Dim>
using Tensor = Eigen::Matrix<double, tensorSize(Dim), 1>;

}  // namespace MaterialPropertyLib
