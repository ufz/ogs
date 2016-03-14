/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "BaseLib/ArrayRef.h"
#include "BaseLib/MatrixRef.h"

namespace MathLib
{
template <unsigned N, typename... FunctionArguments>
using Function = bool (*)(const double t,
                          BaseLib::ArrayRef<const double, N> y,
                          BaseLib::ArrayRef<double, N> ydot,
                          FunctionArguments&... arg);

template <unsigned N, typename... FunctionArguments>
using JacobianFunction = bool (*)(const double t,
                                  BaseLib::ArrayRef<const double, N> y,
                                  BaseLib::ArrayRef<const double, N> ydot,
                                  BaseLib::MatrixRef<double, N, N> jac,
                                  FunctionArguments&... arg);

// This is an internal detail
class FunctionHandles
{
public:
	virtual bool call(const double t, double const* const y,
	                  double* const ydot) = 0;
	virtual bool callJacobian(const double t,
	                          double const* const y,
	                          double const* const ydot,
	                          double* const jac,
	                          BaseLib::StorageOrder order) = 0;

	virtual bool hasJacobian() const = 0;

	virtual unsigned getNumEquations() const = 0;

	virtual ~FunctionHandles() = default;
};
}
