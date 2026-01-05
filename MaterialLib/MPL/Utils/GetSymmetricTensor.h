// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

#include "MaterialLib/MPL/Property.h"
#include "MathLib/KelvinVector.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
constexpr int symmetric_tensor_size =
    MathLib::KelvinVector::kelvin_vector_dimensions(GlobalDim);

template <int GlobalDim>
using SymmetricTensor =
    Eigen::Matrix<double, symmetric_tensor_size<GlobalDim>, 1>;

template <int GlobalDim>
SymmetricTensor<GlobalDim> getSymmetricTensor(
    MaterialPropertyLib::PropertyDataType const& values);
}  // namespace MaterialPropertyLib
