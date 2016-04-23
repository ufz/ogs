/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_ODE_ODESOLVERTYPES_H
#define MATHLIB_ODE_ODESOLVERTYPES_H

#include <functional>
#include <Eigen/Core>

namespace MathLib
{
template <int M, int N>
using MappedMatrix = Eigen::Map<Eigen::Matrix<double, M, N, Eigen::ColMajor>>;

template <int M, int N>
using MappedConstMatrix =
    Eigen::Map<const Eigen::Matrix<double, M, N, Eigen::ColMajor>>;

template <int N>
using MappedVector = MappedMatrix<N, 1>;

template <int N>
using MappedConstVector = MappedConstMatrix<N, 1>;

template <unsigned N>
using Function = std::function<bool(
    const double t, MappedConstVector<N> const& y, MappedVector<N>& ydot)>;

template <unsigned N>
using JacobianFunction = std::function<bool(const double t,
                                            MappedConstVector<N> const& y,
                                            MappedConstVector<N> const& ydot,
                                            MappedMatrix<N, N>& jac)>;

}  // namespace MathLib

#endif  // MATHLIB_ODE_ODESOLVERTYPES_H
