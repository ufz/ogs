/**
 * \author Norihiro Watanabe
 * \date   2012-06-25
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_NONLINEAR_NEWTONRAPHSON_H_
#define MATHLIB_NONLINEAR_NEWTONRAPHSON_H_

#include "LinAlg/VectorNorms.h"

namespace MathLib
{

namespace Nonlinear
{

/**
 * \brief Newton-Raphson method
 */
class NewtonRaphson
{
public:
    /// Default constructor with norm type INFINITY_N, relative tolerance 1e-6
    /// for residual, and the maximum number of iterations 25
    NewtonRaphson();

    /// set a vector norm type
    void setNormType(VecNormType normType) {_normType = normType;}

    /// set the maximum number of iterations
    void setMaxIterations(std::size_t max_itr) {_max_itr = max_itr;}

    /// set the absolute residual tolerance used by the stopping criteria
    void setAbsResidualTolerance(double abs_tol) {_r_abs_tol = abs_tol;}

    /// set the relative residual tolerance used by the stopping criteria
    void setRelResidualTolerance(double rel_tol) {_r_rel_tol = rel_tol;}

    /// set the relative solution increment tolerance used by the stopping criteria
    void setRelDxTolerance(double rel_tol) {_dx_rel_tol = rel_tol;}

    /// print errors during iterations
    void printErrors(bool flag) {_printErrors = flag;}

    /**
     * solve a nonlinear problem
     *
     * \tparam F_RESIDUAL       Function object returning a residual of an equation.
     * The object should have an operator ()(\f$x_k\f$, \f$r_k}\f$).
     * \f$x_k\f$ is current solution. \f$r_k\f$ is calculated residual and an output
     * of this function.
     * \tparam F_DX             Function object returning a solution increment.
     * The object should have an operator ()(\f$x_k\f$, \f$r_k, \f$\Delta x_{k+1}}\f$).
     * \f$\Delta x_{k+1}}\f$ is a calculated solution increment and an output of
     * this function.
     * \tparam T_VALUE          Data type of \f$x_k\f$
     * Both scalar and vector types are available as far as the following conditions
     * are satisfied
     * - T_VALUE has a copy constructor (for non-primitive data types)
     * - MathLib::norm(T_VALUE) exists
     * \param f_residual        Residual function object \f$r(x)\f$
     * \param f_dx              Solution increment function object \f$dx(x)=J^{-1}r\f$
     * \param x0                Initial guess
     * \param x_new             Solution
     * \return true if converged
     */
    template<class F_RESIDUAL, class F_DX, class T_VALUE>
    bool solve(F_RESIDUAL &f_residual, F_DX &f_dx, const T_VALUE &x0, T_VALUE &x_new);

    /// return the number of iterations
    std::size_t getNIterations() const {return _n_iterations; }

    /// return absolute error in the last iteration
    double getAbsResidualError() const {return _r_abs_error; }

    /// return relative error in the last iteration
    double getRelResidualError() const {return _r_rel_error; }

    /// return relative error in the last iteration
    double getRelDxError() const {return _dx_rel_error; }

private:
#ifdef NDEBUG
    /// print out for debugging
    template<class T_VALUE>
    inline void printout(std::ostream& os, std::size_t i, T_VALUE& x_new, T_VALUE& r, T_VALUE& dx);
#endif

    /// vector norm type used to evaluate errors
    VecNormType _normType;
    /// absolute tolerance for residual
    double _r_abs_tol;
    /// relative tolerance for residual
    double _r_rel_tol;
    /// relative tolerance for dx
    double _dx_rel_tol;
    /// the maximum allowed iteration number
    std::size_t _max_itr;
    /// print iteration errors or not
    bool _printErrors;
    /// the number of iterations in the last calculation
    std::size_t _n_iterations;
    /// absolute residual error in the last calculation
    double _r_abs_error;
    /// relative residual error in the last calculation
    double _r_rel_error;
    /// relative dx error in the last calculation
    double _dx_rel_error;
};

} // Nonlinear
} // MathLib

#include "NewtonRaphson-impl.h"

#endif // MATHLIB_NONLINEAR_NEWTONRAPHSON_H_
