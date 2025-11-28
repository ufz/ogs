/*
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 31, 2019, 11:28 AM
 */

#include "FormEigenTensor.h"

#include "MaterialLib/MPL/PropertyType.h"
#include "MathLib/FormattingUtils.h"

namespace MaterialPropertyLib
{

template Eigen::Matrix<double, 1, 1> formEigenTensor<1>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 2, 2> formEigenTensor<2>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 3, 3> formEigenTensor<3>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 4, 4> formEigenTensor<4>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 6, 6> formEigenTensor<6>(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
