/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_ODE_HANDLES_H
#define MATHLIB_ODE_HANDLES_H

#include "OdeSolverTypes.h"

namespace MathLib
{
namespace detail
{
/// Function handles for N equations.
template <unsigned N>
struct Handles : public MathLib::FunctionHandles
{
	using Function = MathLib::Function<N>;
	using JacobianFunction = MathLib::JacobianFunction<N>;

	Handles(Function& f, JacobianFunction& df) : f(f), df(df) {}
	bool call(const double t, const double* const y,
	          double* const ydot) override
	{
		if (f) return f(t, MappedConstVector<N>{y}, MappedVector<N>{ydot});
		return false;
	}

	bool callJacobian(const double t, const double* const y, double* const ydot,
	                  double* const jac) override
	{
		if (df)
			return df(t,
			          MappedConstVector<N>{y},
			          MappedConstVector<N>{ydot},
			          MappedMatrix<N, N>{jac});
		return false;
	}

	bool hasJacobian() const override { return df != nullptr; }
	unsigned getNumEquations() const override { return N; }
	Function f;
	JacobianFunction df;
};

}  // namespace detail
}  // namespace MathLib

#endif  // MATHLIB_ODE_HANDLES_H
