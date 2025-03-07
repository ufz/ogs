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
#include "MathLib/LinAlg/LinearSolverBehaviour.h"
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
    virtual void assemble(std::vector<GlobalVector*> const& x,
                          std::vector<GlobalVector*> const& x_prev,
                          int const process_id) = 0;

    /// \return The global indices for the entries of the global residuum
    /// vector that do not need initial non-equilibrium compensation.
    virtual std::vector<GlobalIndexType>
    getIndicesOfResiduumWithoutInitialCompensation() const = 0;

    /*! Writes the residual at point \c x to \c res.
     *
     * \pre assemble() must have been called before with the same argument \c x.
     *
     * \todo Remove argument \c x.
     */
    virtual void getResidual(GlobalVector const& x,
                             GlobalVector const& x_prev,
                             GlobalVector& res) const = 0;

    /*! Writes the Jacobian of the residual to \c Jac.
     *
     * \pre assemble() must have been called before.
     */
    virtual void getJacobian(GlobalMatrix& Jac) const = 0;

    //! Pre-compute known solutions and possibly store them internally.
    virtual void computeKnownSolutions(GlobalVector const& x,
                                       int const process_id) = 0;

    //! Apply known solutions to the solution vector \c x.
    //! \pre computeKnownSolutions() must have been called before.
    virtual void applyKnownSolutions(GlobalVector& x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ \mathit{Jac} \cdot (-\Delta x) = \mathit{res} \f$.
    //! \pre computeKnownSolutions() must have been called before.
    virtual void applyKnownSolutionsNewton(
        GlobalMatrix& Jac, GlobalVector& res, GlobalVector const& x,
        GlobalVector& minus_delta_x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ \mathit{Jac} \cdot (-\Delta x) = \mathit{res} \f$.
    //! \pre computeKnownSolutions() must have been called before.
    virtual void applyKnownSolutionsPETScSNES(GlobalMatrix& Jac,
                                              GlobalVector& res,
                                              GlobalVector& x) const = 0;

    virtual void updateConstraints(GlobalVector& /*lower*/,
                                   GlobalVector& /*upper*/,
                                   int /*process_id*/) = 0;
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
    virtual void assemble(std::vector<GlobalVector*> const& x,
                          std::vector<GlobalVector*> const& x_prev,
                          int const process_id) = 0;

    /// \return The global indices for the entries of the global residuum
    /// vector that do not need initial non-equilibrium compensation.
    virtual std::vector<GlobalIndexType>
    getIndicesOfResiduumWithoutInitialCompensation() const = 0;

    //! Writes the linearized equation system matrix to \c A.
    //! \pre assemble() must have been called before.
    virtual void getA(GlobalMatrix& A) const = 0;

    //! Writes the linearized equation system right-hand side to \c rhs.
    //! \pre assemble() must have been called before.
    virtual void getRhs(GlobalVector const& x_prev,
                        GlobalVector& rhs) const = 0;

    //! Writes the A_transposed times A into \c A
    //! and also writes A_transposed times rhs into \c rhs
    //! \pre getA() and getRhs must have been called before.
    virtual void getAandRhsNormalized(GlobalMatrix& A,
                                      GlobalVector& rhs) const = 0;

    //! Pre-compute known solutions and possibly store them internally.
    virtual void computeKnownSolutions(GlobalVector const& x,
                                       int const process_id) = 0;

    //! Apply known solutions to the solution vector \c x.
    //! \pre computeKnownSolutions() must have been called before.
    virtual void applyKnownSolutions(GlobalVector& x) const = 0;

    //! Apply known solutions to the linearized equation system
    //! \f$ A \cdot x = \mathit{rhs} \f$.
    //! \pre computeKnownSolutions() must have been called before.
    virtual void applyKnownSolutionsPicard(GlobalMatrix& A, GlobalVector& rhs,
                                           GlobalVector& x) const = 0;

    //! Returns whether the assembled matrix \f$A\f$ has changed and the linear
    //! solver must perform the MathLib::EigenLinearSolver::compute() step.
    virtual MathLib::LinearSolverBehaviour linearSolverNeedsToCompute()
        const = 0;
};

//! @}
}  // namespace NumLib
