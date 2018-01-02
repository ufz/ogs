/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "EquationSystem.h"
#include "Types.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

/*! A System of nonlinear equations.
 *
 * \tparam NLTag  a tag indicating the method used for solving the equation.
 */
template <NonlinearSolverTag NLTag>
class NonlinearSystem;

/*! A System of nonlinear equations to be solved with the Newton-Raphson method.
 *
 * The Newton-Raphson method will iterate the linearized equation
 * \f$ \mathtt{Jac} \cdot (-\Delta x_i) = \mathtt{res} \f$.
 */
template <>
class NonlinearSystem<NonlinearSolverTag::Newton> : public EquationSystem
{
public:
    //! Assembles the linearized equation system at the point \c x.
    //! The linearized system is \f$A(x) \cdot x = b(x)\f$. Here the matrix
    //! \f$A(x)\f$ and the vector \f$b(x)\f$ are assembled.
    virtual void assemble(GlobalVector const& x) = 0;

    /*! Writes the residual at point \c x to \c res.
     *
     * \pre assemble() must have been called before with the same argument \c x.
     *
     * \todo Remove argument \c x.
     */
    virtual void getResidual(GlobalVector const& x,
                             GlobalVector& res) const = 0;

    /*! Writes the Jacobian of the residual to \c Jac.
     *
     * \pre assemble() must have been called before.
     */
    virtual void getJacobian(GlobalMatrix& Jac) const = 0;

    //! Apply known solutions to the solution vector \c x.
    virtual void applyKnownSolutions(GlobalVector& x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ \mathit{Jac} \cdot (-\Delta x) = \mathit{res} \f$.
    virtual void applyKnownSolutionsNewton(GlobalMatrix& Jac, GlobalVector& res,
                                           GlobalVector& minus_delta_x) = 0;
};

/*! A System of nonlinear equations to be solved with the Picard fixpoint
 *  iteration method.
 *
 * The Picard method will iterate the linearized equation
 * \f$ \mathtt{A} \cdot x_i = \mathtt{rhs} \f$.
 */
template <>
class NonlinearSystem<NonlinearSolverTag::Picard> : public EquationSystem
{
public:
    //! Assembles the linearized equation system at the point \c x.
    //! The linearized system is \f$J(x) \cdot \Delta x = (x)\f$. Here the
    //! residual vector \f$r(x)\f$ and its Jacobian \f$J(x)\f$ are assembled.
    virtual void assemble(GlobalVector const& x) = 0;

    //! Writes the linearized equation system matrix to \c A.
    //! \pre assemble() must have been called before.
    virtual void getA(GlobalMatrix& A) const = 0;

    //! Writes the linearized equation system right-hand side to \c rhs.
    //! \pre assemble() must have been called before.
    virtual void getRhs(GlobalVector& rhs) const = 0;

    //! Apply known solutions to the solution vector \c x.
    virtual void applyKnownSolutions(GlobalVector& x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ A \cdot x = \mathit{rhs} \f$.
    virtual void applyKnownSolutionsPicard(GlobalMatrix& A, GlobalVector& rhs,
                                           GlobalVector& x) = 0;
};

//! @}
}
