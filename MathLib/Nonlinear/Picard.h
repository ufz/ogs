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

#ifndef MATHLIB_NONLINEAR_PICARD_H_
#define MATHLIB_NONLINEAR_PICARD_H_

#include "LinAlg/VectorNorms.h"

namespace MathLib
{

namespace Nonlinear
{

/**
 * \brief Picard method (Fixed-point iteration method)
 */
class Picard
{
public:
    /// Default constructor with norm type INFINITY_N, relative tolerance 1e-6,
    /// and the maximum number of iterations 25
    Picard();

    /// set a vector norm type
    void setNormType(VecNormType normType) {_normType = normType;}

    /// set the maximum number of iterations
    void setMaxIterations(std::size_t max_itr) {_max_itr = max_itr;}

    /// set the absolute tolerance used by the stopping criteria
    void setAbsTolerance(double abs_tol) {_abs_tol = abs_tol;}

    /// set the relative tolerance used by the stopping criteria
    void setRelTolerance(double rel_tol) {_rel_tol = rel_tol;}

    /// print errors during iterations
    void printErrors(bool flag) {_printErrors = flag;}

    /**
     * solve a nonlinear problem
     *
     * \param g                 Fixed point function object (x=g(x))
     * \param x0                Initial guess
     * \param x_new             Solution
     * \return true if converged
     */
    template<class F_PROBLEM, class T_VALUE>
    bool solve(F_PROBLEM &g,  const T_VALUE &x0, T_VALUE &x_new);

    /// return the number of iterations
    std::size_t getNIterations() const {return _n_iterations; }

    /// return absolute error in the last iteration
    double getAbsError() const {return _abs_error; }

    /// return relative error in the last iteration
    double getRelError() const {return _rel_error; }

private:
#ifdef NDEBUG
    /// print out for debugging
    template<class T_VALUE>
    inline void printout(std::ostream& os, std::size_t i, T_VALUE& x_new, T_VALUE& dx);
#endif

    /// vector norm type used to evaluate errors
    VecNormType _normType;
    /// absolute tolerance for dx
    double _abs_tol;
    /// relative tolerance for dx
    double _rel_tol;
    /// the maximum allowed iteration number
    std::size_t _max_itr;
    /// print iteration errors or not
    bool _printErrors;
    /// the number of iterations in the last calculation
    std::size_t _n_iterations;
    /// absolute error in the last calculation
    double _abs_error;
    /// relative error in the last calculation
    double _rel_error;
};

} // Nonlinear
} // MathLib

#include "Picard-impl.h"

#endif // MATHLIB_NONLINEAR_PICARD_H_
