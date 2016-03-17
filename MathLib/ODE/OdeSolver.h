/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_ODESOLVER_H
#define MATHLIB_ODESOLVER_H

#include <array>

#include "declarations.h"

namespace MathLib
{
/**
 * ODE solver Interface.
 *
 * This class provides an abstract interface for ODE solvers.
 * It provides type-safe and array-bounds checked access to external
 * ODE solver libraries. However, it is agnostic to the specific solver used.
 */
template <unsigned NumEquations, typename... FunctionArguments>
class OdeSolver
{
public:
	using Arr = std::array<double, NumEquations>;
	using ConstArrRef =
	    Eigen::Map<const Eigen::Matrix<double, NumEquations, 1>>;
	using Function = MathLib::Function<NumEquations, FunctionArguments...>;
	using JacobianFunction =
	    MathLib::JacobianFunction<NumEquations, FunctionArguments...>;

	virtual void setFunction(Function f, JacobianFunction df,
	                         FunctionArguments&... args) = 0;

	virtual void setTolerance(const Arr& abstol, const double reltol) = 0;
	virtual void setTolerance(const double abstol, const double reltol) = 0;

	virtual void setIC(const double t0, const Arr& y0) = 0;

	virtual void preSolve() = 0;
	virtual void solve(const double t) = 0;

	virtual unsigned getNumEquations() const { return NumEquations; }
	virtual ConstArrRef getSolution() const = 0;
	virtual double getTime() const = 0;
	virtual Arr getYDot(const double t, const Arr& y) const = 0;

	virtual ~OdeSolver() = default;
};

}  // namespace MathLib

#endif  // MATHLIB_ODESOLVER_H
