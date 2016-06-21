/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_NONLINEARSYSTEM_H
#define NUMLIB_NONLINEARSYSTEM_H

#include "EquationSystem.h"
#include "Types.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

/*! A System of nonlinear equations.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the
 *                equation.
 * \tparam Vector the type of the solution vector of the equation.
 * \tparam NLTag  a tag indicating the method used for solving the equation.
 */
template <typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class NonlinearSystem;

/*! A System of nonlinear equations to be solved with the Newton-Raphson method.
 *
 * The Newton-Raphson method will iterate the linearized equation
 * \f$ \mathtt{Jac} \cdot (-\Delta x_i) = \mathtt{res} \f$.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the
 *                equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template <typename Matrix, typename Vector>
class NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton>
    : public EquationSystem<Vector>
{
public:
    //! Assembles the residual at the point \c x.
    virtual void assembleResidualNewton(Vector const& x) = 0;

    //! Assembles the Jacobian of the residual at the point \c x.
    virtual void assembleJacobian(Vector const& x) = 0;

    /*! Writes the residual at point \c x to \c res.
     *
     * \pre assembleResidualNewton() must have been called before
     *      with the same argument \c x.
     *
     * \todo Remove argument \c x.
     */
    virtual void getResidual(Vector const& x, Vector& res) const = 0;

    /*! Writes the Jacobian of the residual to \c Jac.
     *
     * \pre assembleJacobian() must have been called before.
     */
    virtual void getJacobian(Matrix& Jac) const = 0;

    //! Apply known solutions to the solution vector \c x.
    virtual void applyKnownSolutions(Vector& x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ \mathit{Jac} \cdot (-\Delta x) = \mathit{res} \f$.
    virtual void applyKnownSolutionsNewton(Matrix& Jac, Vector& res,
                                           Vector& minus_delta_x) = 0;
};

/*! A System of nonlinear equations to be solved with the Picard fixpoint
 *  iteration method.
 *
 * The Picard method will iterate the linearized equation
 * \f$ \mathtt{A} \cdot x_i = \mathtt{rhs} \f$.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the
 *                equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template <typename Matrix, typename Vector>
class NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard>
    : public EquationSystem<Vector>
{
public:
    //! Assembles the linearized equation at point \c x.
    virtual void assembleMatricesPicard(Vector const& x) = 0;

    //! Writes the linearized equation system matrix to \c A.
    virtual void getA(Matrix& A) const = 0;

    //! Writes the linearized equation system right-hand side to \c rhs.
    virtual void getRhs(Vector& rhs) const = 0;

    //! Apply known solutions to the solution vector \c x.
    virtual void applyKnownSolutions(Vector& x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ A \cdot x = \mathit{rhs} \f$.
    virtual void applyKnownSolutionsPicard(Matrix& A, Vector& rhs,
                                           Vector& x) = 0;
};

//! @}
}

#endif  // NUMLIB_NONLINEARSYSTEM_H
