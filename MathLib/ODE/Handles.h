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

namespace MathLib
{
namespace detail
{

/// Function handles for N equations and arbitrary arguments.
template <unsigned N, typename... FunctionArguments>
struct Handles;

/// Function handles for N equations and single argument.
template <unsigned N, typename FunctionArgument>
struct Handles<N, FunctionArgument> : public MathLib::FunctionHandles
{
	using Function = MathLib::Function<N, FunctionArgument>;
	using JacobianFunction = MathLib::JacobianFunction<N, FunctionArgument>;

	Handles(Function& f, JacobianFunction& df, FunctionArgument& arg)
	    : f(f), df(df), _data(arg)
	{
	}

	bool call(const double t, const double* const y,
	          double* const ydot) override
	{
		// looks like f and df could be any callable object with suitable
		// signature
		// consider omission of data pointer and switch to std::function or
		// alike
		if (f)
			return f(t, MappedConstVector<N>{y}, MappedVector<N>{ydot}, _data);
		return true;
	}

	bool callJacobian(const double t, const double* const y, double* const ydot,
	                  double* const jac) override
	{
		if (df)
			return df(t,
			          MappedConstVector<N>{y},
			          MappedVector<N>{ydot},
			          MappedMatrix<N, N>{jac},
			          _data);
		return true;
	}

	bool hasJacobian() const override { return df != nullptr; }
	unsigned getNumEquations() const override { return N; }
private:
	Function f = nullptr;
	JacobianFunction df = nullptr;
	FunctionArgument& _data;
};

/// Function handles for N equations and no arguments.
template <unsigned N>
struct Handles<N> : public MathLib::FunctionHandles
{
	using Function = MathLib::Function<N>;
	using JacobianFunction = MathLib::JacobianFunction<N>;

	Handles(Function& f, JacobianFunction& df) : f(f), df(df) {}
	bool call(const double t, const double* const y,
	          double* const ydot) override
	{
		if (f)
		{
			return f(t, MappedConstVector<N>{y}, MappedVector<N>{ydot});
		}
		return true;
	}

	bool callJacobian(const double t, const double* const y, double* const ydot,
	                  double* const jac) override
	{
		if (df)
			return df(t,
			          MappedConstVector<N>{y},
			          MappedVector<N>{ydot},
			          MappedMatrix<N, N>{jac});
		return true;
	}

	bool hasJacobian() const override { return df != nullptr; }
	unsigned getNumEquations() const override { return N; }
	Function f = nullptr;
	JacobianFunction df = nullptr;
};

}  // namespace detail
}  // namespace MathLib

#endif  // MATHLIB_ODE_HANDLES_H
