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

#include <Eigen/Core>

namespace MathLib
{
template <int M, int N>
using MappedMatrix = Eigen::Map<Eigen::Matrix<double, M, N>>;

template <int M, int N>
using MappedConstMatrix = Eigen::Map<const Eigen::Matrix<double, M, N>>;

template <int N>
using MappedVector = MappedMatrix<N, 1>;

template <int N>
using MappedConstVector = MappedConstMatrix<N, 1>;

template <unsigned N, typename... FunctionArguments>
using Function = bool (*)(const double t,
                          MappedConstVector<N> const y,
                          MappedVector<N> ydot,
                          FunctionArguments&... arg);

template <unsigned N, typename... FunctionArguments>
using JacobianFunction = bool (*)(const double t,
                                  MappedConstVector<N> const y,
                                  MappedVector<N> ydot,
                                  MappedMatrix<N, N> jac,
                                  FunctionArguments&... arg);

// This is an internal detail
class FunctionHandles
{
public:
	virtual bool call(const double t, double const* const y,
	                  double* const ydot) = 0;
	virtual bool callJacobian(const double t,
	                          double const* const y,
	                          double* const ydot,
	                          double* const jac) = 0;

	virtual bool hasJacobian() const = 0;

	virtual unsigned getNumEquations() const = 0;

	virtual ~FunctionHandles() = default;
};
}

#endif  // MATHLIB_ODE_ODESOLVERTYPES_H
