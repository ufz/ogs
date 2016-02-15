#pragma once

#include "Types.h"


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
template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class NonlinearSystem;


/*! A System of nonlinear equations to be solved with the Newton-Raphson method.
 *
 * The Newton-Raphson method will iterate the linearized equation
 * \f$ \mathtt{Jac} \cdot (-\Delta x_i) = \mathtt{res} \f$.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template<typename Matrix, typename Vector>
class NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton>
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
     * \todo make const
     */
    virtual void getResidual(Vector const& x, Vector& res) = 0;

    /*! Writes the Jacobian of the residual to \c Jac.
     *
     * \pre assembleJacobian() must have been called before.
     * \todo make const
     */
    virtual void getJacobian(Matrix& Jac) = 0;

    /*! Check whether this is actually a linear equation system.
     *
     * \remark
     * Depending on its parameters an in general nonlinear equation system
     * can be linear in special cases. With this method it is possible to
     * detect that at runtime and thus save an assembly call.
     */
    virtual bool isLinear() const = 0;

    virtual ~NonlinearSystem() = default;
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
template<typename Matrix, typename Vector>
class NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard>
{
public:
    //! Assembles the linearized eqation at point \c x.
    virtual void assembleMatricesPicard(Vector const& x) = 0;

    //! Writes the linearized equation system matrix to \c A.
    //! \todo make const
    virtual void getA(Matrix& A) = 0;

    //! Writes the linearized equation system right-hand side to \c rhs.
    //! //! \todo make const
    virtual void getRhs(Vector& rhs) = 0;

    /*! Check whether this is actually a linear equation system.
     *
     * \remark
     * Depending on its parameters an in general nonlinear equation system
     * can be linear in special cases. With this method it is possible to
     * detect that at runtime and thus save an assembly call.
     */
    virtual bool isLinear() const = 0;

    virtual ~NonlinearSystem() = default;
};

//! @}

}
