/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_FRACTURE_IDENTITY2_H_
#define MATERIALLIB_FRACTURE_IDENTITY2_H_

#include <Eigen/Eigen>

namespace MaterialLib
{
namespace Fracture
{

template <int DisplacementDim>
struct FractureIdentity2
{
    using VectorType = Eigen::Matrix<double, DisplacementDim, 1>;

    static VectorType const value;
};

}  // namespace Fracture
}  // namespace MaterialLib

#endif // MATERIALLIB_FRACTURE_IDENTITY2_H_
