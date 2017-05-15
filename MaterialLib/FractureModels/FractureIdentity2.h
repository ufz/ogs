/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include "materiallib_export.h"

namespace MaterialLib
{
namespace Fracture
{
/// Identity vector for extracting non-shear components of a fracture stress
/// vector and a relative displacement vector, e.g. \f$\mathbf{m}=(0,1)^T\f$ in
/// 2D case.
template <int DisplacementDim>
struct FractureIdentity2
{
    using VectorType = Eigen::Matrix<double, DisplacementDim, 1>;

    static MATERIALLIB_EXPORT VectorType const value;
};

extern template struct FractureIdentity2<2>;
extern template struct FractureIdentity2<3>;

}  // namespace Fracture
}  // namespace MaterialLib
