/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "ODESolverTypes.h"
#include "FunctionHandles.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
namespace ODE
{
//! \addtogroup ExternalODESolverInterface
//! @{

class CVodeSolverImpl;

/*! ODE solver interfacing with Sundials' CVode.
 *
 * All methods of this class have the same semantics as methods from the
 * ODESolver interface with the same name.
 *
 * There is no implicit information about vector and matrix sizes in this class
 * anymore, rather this class only handles raw pointers. But this is exactly
 * what is necessary when using C libraries like Sundials.
 *
 * \note
 * This class is for internal use only. Therefore all members of this class are
 * protected, namely because ConcreteODESolver derives from this class.
 */
class CVodeSolver
{
protected:
    //! Construct from the given \c config with storage allocated for the given
    //! \c num_equations.
    CVodeSolver(BaseLib::ConfigTree const& config,
                unsigned const num_equations);

    void setTolerance(double const* const abstol, const double reltol);
    void setTolerance(const double abstol, const double reltol);

    void setFunction(std::unique_ptr<detail::FunctionHandles>&& f);

    void setIC(const double t0, double const* const y0);

    void preSolve();
    bool solve(const double t_end);

    double const* getSolution() const;
    double getTime() const;
    void getYDot(const double t,
                 double const* const y,
                 double* const y_dot) const;

    ~CVodeSolver();

private:
    //! pimpl idiom.
    std::unique_ptr<CVodeSolverImpl> _impl;
};

//! @}}

}  // namespace ODE
}  // namespace MathLib
