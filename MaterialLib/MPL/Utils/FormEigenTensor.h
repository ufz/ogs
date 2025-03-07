/*
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 29, 2019, 2:48 PM
 */

#pragma once

#include <Eigen/Core>

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, GlobalDim> formEigenTensor(
    MaterialPropertyLib::PropertyDataType const& values);
}  // namespace MaterialPropertyLib
