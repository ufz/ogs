/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_ODESYSTEM_H
#define NUMLIB_ODESYSTEM_H

#include "NumLib/IndexValueVector.h"
#include "NumLib/DOF/MatrixVectorTraits.h"

#include "Types.h"
#include "EquationSystem.h"

namespace NumLib
{

//! \addtogroup ODESolver
//! @{

/*! ODE system interface.
 *
 * This is the interface an ODE has to implement in order to be solved with this
 * ODE solver library.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 * \tparam ODETag a tag indicating the type of ODE.
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template<typename Matrix, typename Vector, ODESystemTag ODETag, NonlinearSolverTag NLTag>
class ODESystem;


/*! Interface for a first-order implicit quasi-linear ODE.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template<typename Matrix, typename Vector>
class ODESystem<Matrix, Vector,
                ODESystemTag::FirstOrderImplicitQuasilinear,
                NonlinearSolverTag::Picard>
        : public EquationSystem<Vector>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    //! Assemble \c M, \c K and \c b at the provided state (\c t, \c x).
    virtual void assemble(const double t, Vector const& x,
                          Matrix& M, Matrix& K, Vector& b) = 0;

    using Index = typename MathLib::MatrixVectorTraits<Matrix>::Index;

    //! Provides known solutions (Dirichlet boundary conditions) vector for
    //! the ode system at the given time \c t.
    virtual std::vector<NumLib::IndexValueVector<Index>> const*
    getKnownSolutions(double const t) const
    {
        (void)t;
        return nullptr; // by default there are no known solutions
    }
};


/*! Interface for a first-order implicit quasi-linear ODE.
 *
 * ODEs using this interface also provide a Jacobian in addition
 * to the functionality of the Picard-related interface.
 */
template<typename Matrix, typename Vector>
class ODESystem<Matrix, Vector,
                ODESystemTag::FirstOrderImplicitQuasilinear,
                NonlinearSolverTag::Newton>
        : public ODESystem<Matrix, Vector,
                           ODESystemTag::FirstOrderImplicitQuasilinear,
                           NonlinearSolverTag::Picard>
{
public:
    /*! Assemble \f$ \mathtt{Jac} := \partial r/\partial x_N \f$ at the provided state (\c t, \c x).
     *
     * For the meaning of the other parameters refer to the the introductory remarks on
     * \ref concept_time_discretization "time discretization".
     *
     * \remark
     * \parblock
     * The Jacobian will be generally of the following form:
     * \f[ \mathtt{Jac} := \frac{\partial r(x_C, t_C)}{\partial x_N} =
     *  M \cdot \frac{\partial \hat x}{\partial x_N}
     *  + \frac{\partial M}{\partial x_N} \cdot \hat x
     *  + K \cdot \frac{\partial x_C}{\partial x_N}
     *  + \frac{\partial K}{\partial x_N} \cdot x_N
     *  + \frac{\partial b}{\partial x_N},
     *  \f]
     * where \f$ M \f$, \f$ K \f$ and \f$ b \f$ are matrix-valued (vector-valued, respectively)
     * functions that depend on \f$ x_C \f$ and \f$ t_C \f$.
     *
     * Due to the arguments provided to this method its implementation only has to
     * compute the derivatives
     * \f$ \frac{\partial M}{\partial x_N} \cdot \hat x \f$,
     * \f$ \frac{\partial K}{\partial x_N} \cdot x_N    \f$ and
     * \f$ \frac{\partial b}{\partial x_N} \f$.
     * The other terms can be readily taken from the method parameters.
     *
     * In particular for the ForwardEuler time discretization scheme the equation will
     * collapse to
     * \f$ \mathtt{Jac} =
     *  M \cdot \frac{\partial \hat x}{\partial x_N}
     *  \f$
     * since in that scheme \f$ x_N \neq x_C \f$.
     *
     * Of course, the implementation of this method is allowed to compute the Jacobian in a
     * different way, as long as that is consistent with the definition of \f$ \mathtt{Jac} \f$.
     * \endparblock
     */
    virtual void assembleJacobian(const double t, Vector const& x, Vector const& xdot,
                                  const double dxdot_dx, Matrix const& M,
                                  const double dx_dx, Matrix const& K,
                                  Matrix& Jac) = 0;
};

//! @}

}

#endif // NUMLIB_ODESYSTEM_H
