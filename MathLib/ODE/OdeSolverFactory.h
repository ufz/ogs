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
#include "Handles.h"

#include "CVodeSolver.h"

namespace MathLib
{

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
	using Function = typename Interface::Function;
	using JacobianFunction = typename Interface::JacobianFunction;

	void setFunction(Function f, JacobianFunction df,
	                 FunctionArguments&... args) override
	{
		Implementation::setFunction(
		    std::unique_ptr<
		        detail::Handles<NumEquations, FunctionArguments...>>{
		        new detail::Handles<NumEquations, FunctionArguments...>{
		            f, df, args...}});
	}

	void setTolerance(const std::array<double, NumEquations>& abstol,
	                  const double reltol) override
	{
		Implementation::setTolerance(abstol.data(), reltol);
	}

	void setTolerance(const double abstol, const double reltol) override
	{
		Implementation::setTolerance(abstol, reltol);
	}

	void setIC(const double t0,
	           std::array<double, NumEquations> const& y0) override
	{
		Implementation::setIC(t0, y0.data());
	}

	void setIC(const double t0,
	           Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> const& y0) override
	{
		Implementation::setIC(t0, y0.data());
	}

	void preSolve() override { Implementation::preSolve(); }
	void solve(const double t) override { Implementation::solve(t); }
	MappedConstVector<NumEquations> getSolution() const override
	{
		return MappedConstVector<NumEquations>{Implementation::getSolution()};
	}

	double getTime() const override { return Implementation::getTime(); }
	Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> getYDot(
	    const double t, const MappedConstVector<NumEquations>& y) const override
	{
		Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> y_dot;
		Implementation::getYDot(t, y.data(), y_dot.data());
		return y_dot;
	}

private:
	/// instances of this class shall only be constructed by
	/// the friend function listed below
	ConcreteOdeSolver(BaseLib::ConfigTree const& config)
	    : Implementation{config, NumEquations}
	{
	}

	friend std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...>>
	createOdeSolver<NumEquations, FunctionArguments...>(
	    BaseLib::ConfigTree const& config);
};

template <unsigned NumEquations, typename... FunctionArguments>
std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...>> createOdeSolver(
    BaseLib::ConfigTree const& config)
{
	return std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...>>(
	    new ConcreteOdeSolver<NumEquations, CVodeSolver, FunctionArguments...>(
	        config));
}

}  // namespace MathLib

#endif  // MATHLIB_ODE_ODESOLVERFACTORY_H
