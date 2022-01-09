/*
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 04, 2020, 5:20 PM
 */

#pragma once

#include <Eigen/Dense>

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
