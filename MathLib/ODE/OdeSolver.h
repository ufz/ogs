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

#include "OdeSolverTypes.h"

namespace MathLib
{
/**
 * ODE solver Interface.
 *
 * This class provides an abstract interface for ODE solvers.
 * It provides type-safe and array-bounds checked access to external
 * ODE solver libraries. However, it is agnostic to the specific solver used.
 */
template <unsigned NumEquations>
class OdeSolver
{
public:
	using Function = MathLib::Function<NumEquations>;
	using JacobianFunction =
	    MathLib::JacobianFunction<NumEquations>;

	virtual void setFunction(Function f, JacobianFunction df) = 0;

	virtual void setTolerance(const std::array<double, NumEquations>& abstol,
	                          const double reltol) = 0;
	virtual void setTolerance(const double abstol, const double reltol) = 0;

	virtual void setIC(const double t0,
	                   std::initializer_list<double> const& y0) = 0;
	virtual void setIC(
	    const double t0,
	    Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> const& y0) = 0;

	virtual void preSolve() = 0;
	virtual void solve(const double t) = 0;

	virtual unsigned getNumEquations() const { return NumEquations; }
	virtual MappedConstVector<NumEquations> getSolution() const = 0;
	virtual double getTime() const = 0;
	virtual Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> getYDot(
	    const double t, const MappedConstVector<NumEquations>& y) const = 0;

	virtual ~OdeSolver() = default;
};

}  // namespace MathLib

#endif  // MATHLIB_ODESOLVER_H
