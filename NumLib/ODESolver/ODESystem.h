#pragma once

#include "Types.h"


namespace NumLib
{

//! \addtogroup ODESolver
//! @{

/*! ODE system interface.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 * \tparam ODETag a tag indicating the type of ODE.
 * \tparam NLTag  a tag indicating the method used for solving the equation.
 */
template<typename Matrix, typename Vector, ODESystemTag ODETag, NonlinearSolverTag NLTag>
class ODESystem;


/*! Interface for a first-order implicit quasi-linear ODE.
 *
 * Such ODEs take the form
 * \f$ M(x,t)\cdot \dot x + K(x,t) \cdot x - b(x,t)
 *  =: r(\dot x, x, t) \stackrel{!}{=} 0 \f$.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 */
template<typename Matrix, typename Vector>
class ODESystem<Matrix, Vector,
                ODESystemTag::FirstOrderImplicitQuasilinear,
                NonlinearSolverTag::Picard>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    /*! Check whether this is actually a linear equation system.
     *
     * \remark
     * Depending on its parameters an in general nonlinear ODE
     * can be linear in special cases. With this method it is possible to
     * detect that at runtime and thus save an assembly call.
     */
    virtual bool isLinear() const = 0;

    //! Get the number of equations.
    virtual IndexType getNumEquations() const = 0;

    //! Assemble \c M, \c K and \c b at the state (\c t, \c x).
    virtual void assemble(const double t, Vector const& x,
                          Matrix& M, Matrix& K, Vector& b) = 0;

    virtual ~ODESystem() = default;
};

template<typename Matrix, typename Vector>
class ODESystem<Matrix, Vector,
                ODESystemTag::FirstOrderImplicitQuasilinear,
                NonlinearSolverTag::Newton>
        : public ODESystem<Matrix, Vector,
                           ODESystemTag::FirstOrderImplicitQuasilinear,
                           NonlinearSolverTag::Picard>
{
public:
    /*! Assemble \f$ \mathtt{Jac} := \partial r/\partial x \f$ at the state (\c t, \c x).
     *
     * For the meaning of the other parameters refer to the the introductory remarks on
     * \ref concept_time_discretization "time discretization".
     */
    virtual void assembleJacobian(const double t, Vector const& x, Vector const& xdot,
                                  const double dxdot_dx, const double dx_dx,
                                  Matrix& Jac) = 0;
};

//! @}

}
