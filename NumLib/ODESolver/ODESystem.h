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

#include "Types.h"


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


/*! Interface for a first-order implicit quasi-linear ODE.
 *
 * ODEs using this interface also provide a Jacobian.
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
    /*! Assemble \f$ \mathtt{Jac} := \partial r/\partial x \f$ at the state (\c t, \c x).
     *
     * For the meaning of the other parameters refer to the the introductory remarks on
     * \ref concept_time_discretization "time discretization".
     *
     * TODO document how to assemble the Jacobian!
     */
    virtual void assembleJacobian(const double t, Vector const& x, Vector const& xdot,
                                  const double dxdot_dx, Matrix const& M,
                                  const double dx_dx, Matrix const& K,
                                  Matrix& Jac) = 0;
};

//! @}

}

#endif // NUMLIB_ODESYSTEM_H
