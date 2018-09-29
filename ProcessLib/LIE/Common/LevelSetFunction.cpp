/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "LevelSetFunction.h"

#include <boost/math/special_functions/sign.hpp>

#include "FractureProperty.h"

namespace
{
// Heaviside step function
inline double Heaviside(double v)
{
    return (v < 0.0) ? 0.0 : 1.0;
}

}  // namespace

namespace ProcessLib
{
namespace LIE
{
double calculateLevelSetFunction(FractureProperty const& frac, double const* x_)
{
    Eigen::Map<Eigen::Vector3d const> x(x_, 3);
    return Heaviside(
        boost::math::sign(frac.normal_vector.dot(x - frac.point_on_fracture)));
}

}  // namespace LIE
}  // namespace ProcessLib
