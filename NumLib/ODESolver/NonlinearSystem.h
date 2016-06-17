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

#include "Types.h"
#include "EquationSystem.h"


namespace NumLib
{

//! \addtogroup ODESolver
//! @{

/*! A System of nonlinear equations.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 * \tparam NLTag  a tag indicating the method used for solving the equation.
 */
template<NonlinearSolverTag NLTag>
class NonlinearSystem;


/*! A System of nonlinear equations to be solved with the Newton-Raphson method.
 *
 * The Newton-Raphson method will iterate the linearized equation
 * \f$ \mathtt{Jac} \cdot (-\Delta x_i) = \mathtt{res} \f$.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template <>
class NonlinearSystem<NonlinearSolverTag::Newton>
        : public EquationSystem
{
public:
    //! Assembles the residual at the point \c x.
    virtual void assembleResidualNewton(GlobalVector const& x) = 0;

    //! Assembles the Jacobian of the residual at the point \c x.
    virtual void assembleJacobian(GlobalVector const& x) = 0;

    /*! Writes the residual at point \c x to \c res.
     *
     * \pre assembleResidualNewton() must have been called before
     *      with the same argument \c x.
     *
     * \todo Remove argument \c x.
     */
    virtual void getResidual(GlobalVector const& x, GlobalVector& res) const = 0;

    /*! Writes the Jacobian of the residual to \c Jac.
     *
     * \pre assembleJacobian() must have been called before.
     */
    virtual void getJacobian(GlobalMatrix& Jac) const = 0;

    //! Apply known solutions to the solution vector \c x.
    virtual void applyKnownSolutions(GlobalVector& x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ \mathit{Jac} \cdot (-\Delta x) = \mathit{res} \f$.
    virtual void applyKnownSolutionsNewton(
            GlobalMatrix& Jac, GlobalVector& res, GlobalVector& minus_delta_x) = 0;
};

/*! A System of nonlinear equations to be solved with the Picard fixpoint
 *  iteration method.
 *
 * The Picard method will iterate the linearized equation
 * \f$ \mathtt{A} \cdot x_i = \mathtt{rhs} \f$.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template <>
class NonlinearSystem<NonlinearSolverTag::Picard>
        : public EquationSystem
{
public:
    //! Assembles the linearized equation at point \c x.
    virtual void assembleMatricesPicard(GlobalVector const& x) = 0;

    //! Writes the linearized equation system matrix to \c A.
    virtual void getA(GlobalMatrix& A) const = 0;

    //! Writes the linearized equation system right-hand side to \c rhs.
    virtual void getRhs(GlobalVector& rhs) const = 0;

    //! Apply known solutions to the solution vector \c x.
    virtual void applyKnownSolutions(GlobalVector& x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ A \cdot x = \mathit{rhs} \f$.
    virtual void applyKnownSolutionsPicard(
            GlobalMatrix& A, GlobalVector& rhs, GlobalVector& x) = 0;
};

//! @}

}

#endif // NUMLIB_NONLINEARSYSTEM_H
