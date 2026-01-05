// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <vector>

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
