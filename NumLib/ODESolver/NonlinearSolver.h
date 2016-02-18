/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_NONLINEARSOLVER_H
#define NUMLIB_NONLINEARSOLVER_H

#include <logog/include/logog.hpp>

#include "Types.h"
#include "NonlinearSystem.h"


namespace NumLib
{

//! \addtogroup ODESolver
//! @{

/*! Find a solution to a nonlinear equation.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 * \tparam NLTag  a tag indicating the method used for solving the equation.
 */
template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class NonlinearSolver;


/*! Find a solution to a nonlinear equation using the Newton-Raphson method.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template<typename Matrix, typename Vector>
class NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Newton> final
{
public:
    //! Type of the nonlinear equation system to be solved.
    using System = NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton>;

    /*! Constructs a new instance.
     *
     * \param tol     the tolerance of the solver. \todo Be more specific about that!
     * \param maxiter the maximum number of iterations used to solve the equation.
     */
    explicit
    NonlinearSolver(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    /*! Only assemble the equation system.
     *
     * \note This method is needed to preload CrankNicolson time discretization scheme.
     *       It is not used for the general solver steps; in those only the solve() method
     *       is needed.
     *
     * \param sys the equation system to be assembled.
     * \param x   the state at which the equation system will be assembled.
     */
    void assemble(System& sys, Vector const& x) const;

    /*! Assemble and solve the equation system.
     *
     * \param sys the equation system to be solved.
     * \param x   in: the initial guess, out: the solution.
     *
     * \retval true if the equation system could be solved
     * \retval false otherwise
     */
    bool solve(System& sys, Vector& x);

private:
    const double _tol;       //!< tolerance of the solver
    const unsigned _maxiter; //!< maximum number of iterations

    Vector _res;             //!< The residual.
    Matrix _J;               //!< The Jacobian of the residual.
    Vector _minus_delta_x;   //!< The Newton-Raphson method solves the linearized equation
                             //!< \f$ J \cdot (-\Delta x) = r \f$ repeatedly.
    double const _alpha = 1; //!< Damping factor. \todo Add constructor parameter.
};


/*! Find a solution to a nonlinear equation using the Picard fixpoint iteration method.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template<typename Matrix, typename Vector>
class NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Picard> final
{
public:
    //! Type of the nonlinear equation system to be solved.
    using System = NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard>;

    /*! Constructs a new instance.
     *
     * \param tol     the tolerance of the solver. \todo Be more specific about that!
     * \param maxiter the maximum number of iterations used to solve the equation.
     */
    explicit
    NonlinearSolver(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    /*! Only assemble the equation system.
     *
     * \note This method is needed to preload CrankNicolson time discretization scheme.
     *       It is not used for the general solver steps; in those only the solve() method
     *       is needed.
     *
     * \param sys the equation system to be assembled.
     * \param x   the state at which the equation system will be assembled.
     */
    void assemble(System& sys, Vector const& x) const;

    /*! Assemble and solve the equation system.
     *
     * \param sys the equation system to be solved.
     * \param x   in: the initial guess, out: the solution.
     *
     * \retval true if the equation system could be solved
     * \retval false otherwise
     */
    bool solve(System& sys, Vector& x);

private:
    const double _tol;       //!< tolerance of the solver
    const unsigned _maxiter; //!< maximum number of iterations

    Matrix _A;     //!< \c Matrix describing the linearized system.
    Vector _rhs;   //!< \c Vector describing the linearized system.
    Vector _x_new; //!< \c Vector to store solutions of \f$ A \cdot x = \mathit{rhs} \f$.
};

//! @}

}

#include "NonlinearSolver-impl.h"

#endif // NUMLIB_NONLINEARSOLVER_H
