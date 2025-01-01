/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "EquationSystem.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumLib/IndexValueVector.h"
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

    /// \return The global indices for the entries of the global residuum
    /// vector that do not need initial non-equilibrium compensation.
    virtual std::vector<GlobalIndexType>
    getIndicesOfResiduumWithoutInitialCompensation() const
    {
        return {};
    }

    //! Assemble \c M, \c K and \c b at the provided state (\c t, \c x).
    virtual void assemble(const double t, double const dt,
                          std::vector<GlobalVector*> const& x,
                          std::vector<GlobalVector*> const& x_prev,
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

    virtual void updateConstraints(GlobalVector& /*lower*/,
                                   GlobalVector& /*upper*/,
                                   int const /*process_id*/){};

    //! Indicates whether the assembled matrices change only upon timestep
    //! change. This enables some optimizations in the linear solver, but the
    //! user must be 100 % sure that she uses this option correctly.
    virtual bool shouldLinearSolverComputeOnlyUponTimestepChange() const
    {
        return false;
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

    /*! Assemble \c b and the Jacobian
     * \f$ \mathtt{Jac} := \partial r/\partial x_N \f$
     * at the provided state (\c t, \c x).
     *
     * For the meaning of the other parameters refer to the the introductory
     * remarks on
     * \ref concept_time_discretization "time discretization".
     */
    virtual void assembleWithJacobian(const double t, double const dt,
                                      std::vector<GlobalVector*> const& x,
                                      std::vector<GlobalVector*> const& x_prev,
                                      int const process_id, GlobalVector& b,
                                      GlobalMatrix& Jac) = 0;
};

//! @}
}  // namespace NumLib
