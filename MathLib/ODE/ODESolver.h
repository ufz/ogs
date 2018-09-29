/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>

#include "ODESolverTypes.h"

namespace MathLib
{
namespace ODE
{
//! \addtogroup ExternalODESolverInterface
//! @{

/**
 * ODE solver Interface.
 *
 * This class provides an abstract interface for ODE solvers.
 * It provides type-safe and array-bounds checked access to external
 * ODE solver libraries. However, it is agnostic to the specific solver used.
 *
 * The ODEs solved using this class are first-order ODEs that are given in
 * explicit form, i.e. \f[ \dot y = f(t, y). \f]
 *
 * \tparam NumEquations number of equations in the ODE system.
 */
template <unsigned NumEquations>
class ODESolver
{
public:
    /*! Sets functions that compute \f$\dot y\f$ and
     *  the Jacobian \f$\partial \dot y/\partial y\f$.
     *
     * If no Jacobian function shall be set, \c nullptr can be passed fo \c df.
     *
     * \remark
     * solve() cannot be directly called after this method, rather preSolve()
     * has to be called first!
     */
    virtual void setFunction(Function<NumEquations> f,
                             JacobianFunction<NumEquations> df) = 0;

    /*! Sets the tolerances for the ODE solver.
     *
     * \param abstol absolute tolerance, one value for all equations.
     * \param reltol relative tolerance.
     *
     * \remark
     * solve() cannot be directly called after this method, rather preSolve()
     * has to be called first!
     */
    virtual void setTolerance(const double abstol, const double reltol) = 0;

    /*! Sets the tolerances for the ODE solver.
     *
     * \param abstol absolute tolerance, one value each equation.
     * \param reltol relative tolerance.
     *
     * \remark
     * solve() cannot be directly called after this method, rather preSolve()
     * has to be called first!
     */
    virtual void setTolerance(const std::array<double, NumEquations>& abstol,
                              const double reltol) = 0;

    /*! Sets the conditions.
     *
     * \param t0 initial time.
     * \param y0 initial values.
     *
     * \remark
     * solve() cannot be directly called after this method, rather preSolve()
     * has to be called first!
     */
    virtual void setIC(const double t0,
                       std::initializer_list<double> const& y0) = 0;

    //! \overload
    virtual void setIC(
        const double t0,
        Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> const& y0) = 0;

    /*! Finishes setting up the ODE solver, makes it ready to solve the provided
     * ODE.
     *
     * This method applies settings to the ODE solver, hence it has to be called
     * after calling setters.
     *
     * \note
     * preSolve() has to be called once before calling solve, it is not
     * necessary
     * to call it after each setter.
     */
    virtual void preSolve() = 0;

    /*! Solves the ODE from the set initial condition to time \c t.
     *
     * \returns true or false indicating whether solving succeeded.
     *
     * \pre preSolve() has to be called before this method.
     */
    virtual bool solve(const double t) = 0;

    //! Returns the number of equations.
    virtual unsigned getNumberOfEquations() const { return NumEquations; }
    //! Returns the solution vector \c y
    virtual MappedConstVector<NumEquations> getSolution() const = 0;

    /*! Returns the time that the solver has reached.
     *
     * The return value should be equal to the time \c t passed to solve() if
     * everything went fine.
     */
    virtual double getTime() const = 0;

    /*! Computes \f$ \dot y = f(t,y) \f$.
     *
     * This method is provided for convenience only.
     */
    virtual Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> getYDot(
        const double t, const MappedConstVector<NumEquations>& y) const = 0;

    virtual ~ODESolver() = default;
};

//! @}

}  // namespace ODE
}  // namespace MathLib
