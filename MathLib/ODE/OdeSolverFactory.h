/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_ODE_ODESOLVERFACTORY_H
#define MATHLIB_ODE_ODESOLVERFACTORY_H

#include <memory>

#include "BaseLib/ConfigTree.h"

#include "OdeSolver.h"

#include "CVodeSolver.h"

namespace MathLib
{
namespace detail
{
template <unsigned N, typename... FunctionArguments>
struct Handles;

template <unsigned N, typename FunctionArgument>
struct Handles<N, FunctionArgument> : public MathLib::FunctionHandles
{
	using Function = MathLib::Function<N, FunctionArgument>;
	using JacobianFunction = MathLib::JacobianFunction<N, FunctionArgument>;

	bool call(const double t, const double* const y,
	          double* const ydot) override
	{
		// looks like f and df could be any callable object with suitable
		// signature
		// consider omission of data pointer and switch to std::function or
		// alike
		if (f)
			return f(t,
			         BaseLib::ArrayRef<const double, N>{y},
			         BaseLib::ArrayRef<double, N>{ydot},
			         *_data);
		return true;
	}

	bool callJacobian(const double t, const double* const y,
	                  const double* const ydot, double* const jac,
	                  BaseLib::StorageOrder order) override
	{
		if (df)
			return df(t,
			          BaseLib::ArrayRef<const double, N>{y},
			          BaseLib::ArrayRef<const double, N>{ydot},
			          BaseLib::MatrixRef<double, N, N>{jac, order},
			          *_data);
		return true;
	}

	bool hasJacobian() const override { return df != nullptr; }
	unsigned getNumEquations() const override { return N; }
	void setArguments(FunctionArgument* arg)
	{
		assert(arg != nullptr);
		_data = arg;
	}

	// TODO: make private
	Function f = nullptr;
	JacobianFunction df = nullptr;

private:
	FunctionArgument* _data = nullptr;
};

template <unsigned N>
struct Handles<N> : public MathLib::FunctionHandles
{
	using Function = MathLib::Function<N>;
	using JacobianFunction = MathLib::JacobianFunction<N>;

	bool call(const double t, const double* const y,
	          double* const ydot) override
	{
		if (f)
			return f(t,
			         BaseLib::ArrayRef<const double, N>{y},
			         BaseLib::ArrayRef<double, N>{ydot});
		return true;
	}

	bool callJacobian(const double t, const double* const y,
	                  const double* const ydot, double* const jac,
	                  BaseLib::StorageOrder order) override
	{
		if (df)
			return df(t,
			          BaseLib::ArrayRef<const double, N>{y},
			          BaseLib::ArrayRef<const double, N>{ydot},
			          BaseLib::MatrixRef<double, N, N>{jac, order});
		return true;
	}

	bool hasJacobian() const override { return df != nullptr; }
	unsigned getNumEquations() const override { return N; }
	void setArguments() const {}
	Function f = nullptr;
	JacobianFunction df = nullptr;
};

}  // namespace detail

template <unsigned NumEquations, typename... FunctionArguments>
std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...>> createOdeSolver(
    BaseLib::ConfigTree const& config);

/**
 * ODE solver with a bounds-safe interface.
 *
 * This class makes contact between the abstract \c OdeSolver interface and a
 * certain solver \c Implementation.
 *
 * The interface of this class inherits the array bounds checking from \c
 * OdeSolver.
 * Its methods forward calls to the \c Implementation erasing array bounds info
 * by
 * passing \c std::array as raw pointer.
 *
 * This way the \c Implementation does not need to be templated.
 */
template <unsigned NumEquations, typename Implementation,
          typename... FunctionArguments>
class ConcreteOdeSolver final
    : public OdeSolver<NumEquations, FunctionArguments...>,
      private Implementation
{
public:
	using Interface = OdeSolver<NumEquations, FunctionArguments...>;
	using Arr = typename Interface::Arr;
	using ConstArrRef = typename Interface::ConstArrRef;
	using Function = typename Interface::Function;
	using JacobianFunction = typename Interface::JacobianFunction;

	void init() override
	{
		Implementation::init(NumEquations);
		Implementation::setFunction(&_handles);
	}

	void setTolerance(const Arr& abstol, const double reltol) override
	{
		Implementation::setTolerance(abstol.data(), reltol);
	}

	void setTolerance(const double abstol, const double reltol) override
	{
		Implementation::setTolerance(abstol, reltol);
	}

	void setFunction(Function f, JacobianFunction df,
	                 FunctionArguments*... args) override
	{
		_handles.f = f;
		_handles.df = df;
		_handles.setArguments(args...);
	}

	void setIC(const double t0, const Arr& y0) override
	{
		Implementation::setIC(t0, y0.data());
	}

	void preSolve() { Implementation::preSolve(); }
	void solve(const double t) override { Implementation::solve(t); }
	ConstArrRef getSolution() const override
	{
		return ConstArrRef(Implementation::getSolution());
	}

	double getTime() const override { return Implementation::getTime(); }
	Arr getYDot(const double t, const Arr& y) const override
	{
		Arr ydot;
		Implementation::getYDot(t, y.data(), ydot.data());
		return ydot;
	}

private:
	/// instances of this class shall only be constructed by
	/// the friend function listed below
	ConcreteOdeSolver(BaseLib::ConfigTree const& config)
	    : Implementation{config}
	{
	}

	detail::Handles<NumEquations, FunctionArguments...> _handles;

	friend std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...>>
	createOdeSolver<NumEquations, FunctionArguments...>(
	    BaseLib::ConfigTree const& config);
};

template <unsigned NumEquations, typename... FunctionArguments>
std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...>> createOdeSolver(
    BaseLib::ConfigTree const& config)
{
	return std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...>>(
	    new ConcreteOdeSolver<NumEquations, CVodeSolverInternal,
	                          FunctionArguments...>(config));
}

}  // namespace MathLib

#endif  // MATHLIB_ODE_ODESOLVERFACTORY_H
