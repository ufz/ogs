/*
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 29, 2019, 2:48 PM
 */

#pragma once

#include <Eigen/Dense>

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, 1> formEigenVector(
    MaterialPropertyLib::PropertyDataType const& values);
}  // namespace MaterialPropertyLib
