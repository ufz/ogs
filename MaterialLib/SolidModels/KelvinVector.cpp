/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "KelvinVector.h"

namespace MaterialLib
{
namespace SolidModels
{
template <>
double Invariants<6>::determinant(Eigen::Matrix<double, 6, 1> const& v)
{
    return v(0) * v(1) * v(2) + v(3) * v(4) * v(5) / std::sqrt(2.) -
           v(3) * v(3) * v(2) / 2. - v(4) * v(4) * v(0) / 2. -
           v(5) * v(5) * v(1) / 2.;
}

template <>
double Invariants<4>::determinant(Eigen::Matrix<double, 4, 1> const& v)
{
    return v(2) * (v(0) * v(1) - v(3) * v(3) / 2.);
}

template <>
Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> inverse(
    Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> const& v)
{
    assert(Invariants<4>::determinant(v) != 0);

    Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> inv;
    inv(0) = v(1) * v(2);
    inv(1) = v(0) * v(2);
    inv(2) = v(0) * v(1) - v(3) * v(3) / 2.;
    inv(3) = -v(3) * v(2);
    return inv / Invariants<4>::determinant(v);
}

template <>
Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> inverse(
    Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> const& v)
{
    assert(Invariants<6>::determinant(v) != 0);

    Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> inv;
    inv(0) = v(1) * v(2) - v(4) * v(4) / 2.;
    inv(1) = v(0) * v(2) - v(5) * v(5) / 2.;
    inv(2) = v(0) * v(1) - v(3) * v(3) / 2.;
    inv(3) = v(4) * v(5) / std::sqrt(2.) - v(3) * v(2);
    inv(4) = v(3) * v(5) / std::sqrt(2.) - v(4) * v(0);
    inv(5) = v(4) * v(3) / std::sqrt(2.) - v(1) * v(5);
    return inv / Invariants<6>::determinant(v);
}

}  // namespace SolidModels
}  // namespace MaterialLib
