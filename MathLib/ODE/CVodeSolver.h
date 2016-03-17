/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_CVODESOLVER_H
#define MATHLIB_CVODESOLVER_H

#include "BaseLib/ConfigTree.h"

#include "declarations.h"

namespace MathLib
{
class CVodeSolverImpl;

/**
 * ODE solver, general, pointer based implementation. No implicit bounds
 * checking
 *
 * For internal use only.
 */
class CVodeSolverInternal
{
protected:
	CVodeSolverInternal(BaseLib::ConfigTree const& config,
	                    unsigned const num_equations);

	void setTolerance(double const* const abstol, const double reltol);
	void setTolerance(const double abstol, const double reltol);

	void setFunction(std::unique_ptr<FunctionHandles>&& f);

	void setIC(const double t0, double const* const y0);

	void preSolve();
	void solve(const double t_end);

	double const* getSolution() const;
	double getTime() const;
	void getYDot(const double t,
	             double const* const y,
	             double* const y_dot) const;

	~CVodeSolverInternal();

private:
	std::unique_ptr<CVodeSolverImpl>
	    _impl;  ///< pimpl idiom hides sundials headers.
};

}  // namespace MathLib

#endif  // MATHLIB_CVODESOLVER_H
