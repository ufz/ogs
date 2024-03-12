/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VectorizedTensor.h"

namespace MathLib::VectorizedTensor
{

bool isTensorConvertibleTo1d(Eigen::Matrix3d const& tensor)
{
    return (tensor(0, 1) == 0. && tensor(1, 0) == 0. &&  //
            tensor(0, 2) == 0. && tensor(2, 0) == 0. &&  //
            tensor(1, 2) == 0. && tensor(2, 1) == 0.);
}

bool isTensorConvertibleTo2d(Eigen::Matrix3d const& tensor)
{
    return (tensor(0, 2) == 0. && tensor(1, 2) == 0. &&  //
            tensor(2, 0) == 0. && tensor(2, 1) == 0.);
}

}  // namespace MathLib::VectorizedTensor
