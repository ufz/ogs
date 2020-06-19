/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumLib/IndexValueVector.h"

#include "EquationSystem.h"
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
 * \tparam ODETag a tag indicating the type of ODE.
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template <ODESystemTag ODETag, NonlinearSolverTag NLTag>
class ODESystem;

/*! Interface for a first-order implicit quasi-linear ODE.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template <>
class ODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                NonlinearSolverTag::Picard> : public EquationSystem
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag =
        ODESystemTag::FirstOrderImplicitQuasilinear;

    //! Calls process' pre-assembly with the provided state (\c t, \c dt, \c x).
    virtual void preAssemble(const double t, double const dt,
                             GlobalVector const& x) = 0;

    //! Assemble \c M, \c K and \c b at the provided state (\c t, \c x).
    virtual void assemble(const double t, double const dt,
                          std::vector<GlobalVector*> const& x,
                          std::vector<GlobalVector*> const& xdot,
                          int const process_id, GlobalMatrix& M,
                          GlobalMatrix& K, GlobalVector& b) = 0;

    using Index = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;

    //! Provides known solutions (Dirichlet boundary conditions) vector for
    //! the ode system at the given time \c t and \c process_id.
    virtual std::vector<NumLib::IndexValueVector<Index>> const*
    getKnownSolutions(double const /*t*/, GlobalVector const& /*x*/,
                      int const /*process_id*/) const
    {
        return nullptr;  // by default there are no known solutions
    }
};

/*! Interface for a first-order implicit quasi-linear ODE.
 *
 * ODEs using this interface also provide a Jacobian in addition
 * to the functionality of the Picard-related interface.
 */
template <>
class ODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                NonlinearSolverTag::Newton>
    : public ODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                       NonlinearSolverTag::Picard>
{
public:
    //! Calls process' pre-assembly with the provided state (\c t, \c dt, \c x).
    void preAssemble(const double t, double const dt,
                     GlobalVector const& x) override = 0;

    /*! Assemble \c M, \c K, \c b and the Jacobian
     * \f$ \mathtt{Jac} := \partial r/\partial x_N \f$
     * at the provided state (\c t, \c x).
     *
     * For the meaning of the other parameters refer to the the introductory
     * remarks on
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
     * where \f$ M \f$, \f$ K \f$ and \f$ b \f$ are matrix-valued
     * (vector-valued, respectively)
     * functions that depend on \f$ x_C \f$ and \f$ t_C \f$.
     *
     * Due to the arguments provided to this method its implementation only has
     * to
     * compute the derivatives
     * \f$ \frac{\partial M}{\partial x_N} \cdot \hat x \f$,
     * \f$ \frac{\partial K}{\partial x_N} \cdot x_N    \f$ and
     * \f$ \frac{\partial b}{\partial x_N} \f$.
     * The other terms can be readily taken from the method parameters.
     *
     * In particular for the ForwardEuler time discretization scheme the
     * equation will
     * collapse to
     * \f$ \mathtt{Jac} =
     *  M \cdot \frac{\partial \hat x}{\partial x_N}
     *  \f$
     * since in that scheme \f$ x_N \neq x_C \f$.
     *
     * Of course, the implementation of this method is allowed to compute the
     * Jacobian in a
     * different way, as long as that is consistent with the definition of \f$
     * \mathtt{Jac} \f$.
     * \endparblock
     */
    virtual void assembleWithJacobian(const double t, double const dt,
                                      std::vector<GlobalVector*> const& x,
                                      std::vector<GlobalVector*> const& xdot,
                                      const double dxdot_dx, const double dx_dx,
                                      int const process_id, GlobalMatrix& M,
                                      GlobalMatrix& K, GlobalVector& b,
                                      GlobalMatrix& Jac) = 0;
};

//! @}
}  // namespace NumLib
