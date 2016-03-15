/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_ODE_DECLARATIONS_H
#define MATHLIB_ODE_DECLARATIONS_H

#include <Eigen/Core>

namespace MathLib
{
enum class StorageOrder
{
	ColumnMajor,
	RowMajor
};

template <unsigned N, typename... FunctionArguments>
using Function = bool (*)(const double t,
                          Eigen::Map<const Eigen::Matrix<double, N, 1>> const y,
                          Eigen::Map<Eigen::Matrix<double, N, 1>> ydot,
                          FunctionArguments&... arg);

template <unsigned N, typename... FunctionArguments>
using JacobianFunction =
    bool (*)(const double t,
             Eigen::Map<const Eigen::Matrix<double, N, 1>> const y,
             Eigen::Map<Eigen::Matrix<double, N, 1>> ydot,
             Eigen::Map<Eigen::Matrix<double, N, N>> jac,
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
	                          double* const jac,
	                          StorageOrder order) = 0;

	virtual bool hasJacobian() const = 0;

	virtual unsigned getNumEquations() const = 0;

	virtual ~FunctionHandles() = default;
};
}

#endif  // MATHLIB_ODE_DECLARATIONS_H
