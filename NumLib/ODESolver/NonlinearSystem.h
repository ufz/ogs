// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "EquationSystem.h"
#include "MathLib/LinAlg/LinAlgEnums.h"
#include "MathLib/LinAlg/LinearSolverBehaviour.h"
#include "NumLib/Exceptions.h"
#include "Types.h"

namespace NumLib
{
/*! Runs \c assemble and reconciles a NumLib::AssemblyException thrown by it
 * across all MPI ranks.
 *
 * \param assemble a callable that runs the assembly.
 * \param failure_message logged when the assembly failed.
 * \param log_level level at which the failure is logged; the default is the
 * level of the Newton solver's fatal assembly abort.
 *
 * \return Whether the assembly succeeded on every rank.
 *
 * Do not remove the reconciliation across ranks: a rank that carried on while
 * another one aborted would block in the collective calls of the subsequent
 * linear solve and deadlock the simulation. The reconciliation itself is not
 * covered by a test: the unit tests in
 * Tests/NumLib/TestNonlinearSolverAssemblyException.cpp run on a single rank,
 * and they are disabled in PETSc builds altogether (issue #1989).
 */
template <typename Assemble>
bool assembledOnAllRanks(
    Assemble const& assemble,
    std::string_view const& failure_message,
    spdlog::level::level_enum log_level = spdlog::level::err)
{
    bool mpi_rank_assembly_ok = true;
    try
    {
        assemble();
    }
    catch (AssemblyException const& e)
    {
        BaseLib::console->log(log_level, "{:s} Reason: {:s}", failure_message,
                              e.what());
        mpi_rank_assembly_ok = false;
    }
    return !BaseLib::MPI::anyOf(!mpi_rank_assembly_ok);
}
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

    virtual void setReleaseNodalForces(GlobalVector const* /*r_neq*/,
                                       int const /*process_id*/) = 0;

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

    //! Apply known solutions to the residual \c res and Jacobian \c Jac at
    //! the iterate \c x passed to PETSc SNES: the Dirichlet rows of \c Jac
    //! become identity rows and \f$ \mathit{res}_i = x_i - g_i \f$ at the
    //! Dirichlet DOFs. The latter vanishes if \c x fulfils the Dirichlet BCs,
    //! but keeps a deviation from them visible to SNES, e.g., one caused by
    //! a trust region step.
    //!
    //! \c x is not modified, \c scratch is overwritten at the Dirichlet DOFs.
    //!
    //! \pre computeKnownSolutions() must have been called before, and \c x
    //! must be locally accessible.
    virtual void applyKnownSolutionsPETScSNES(GlobalMatrix& Jac,
                                              GlobalVector& res,
                                              GlobalVector const& x,
                                              GlobalVector& scratch) const = 0;

    //! Writes \f$ P\,(\mathit{res} - J V) \f$ to \c w. \f$ J \f$ is the
    //! Jacobian assembled by the last assemble() call,
    //! \f$ V_k = x_k - g_k \f$ at the known solutions' DOFs \f$ k \f$ and zero
    //! elsewhere, and \f$ P \f$ zeroes the entries at the known solutions'
    //! DOFs.
    //!
    //! For \c res the residual at \c x, this is the residual linearized at
    //! the state fulfilling the Dirichlet BCs, restricted to the free DOFs:
    //! it contains neither reaction forces nor the Dirichlet rows of the
    //! linear system, and it is the same for all linear algebra backends.
    //!
    //! \pre assemble() must have been called at \c x, computeKnownSolutions()
    //! before, and \c x must be locally accessible. \c w must have the layout
    //! of \c res.
    //!
    //! \return \f$ \|V\|_\infty \f$, the largest deviation of \c x from the
    //! Dirichlet values on all MPI ranks.
    virtual double getDirichletFreeResidual(GlobalVector const& x,
                                            GlobalVector const& res,
                                            GlobalVector& w) = 0;

    //! Sets the entries of \c v at the known solutions' DOFs to zero, i.e.,
    //! applies the projection \f$ P \f$ to the free DOFs.
    //!
    //! \pre computeKnownSolutions() must have been called before.
    virtual void projectToFreeDofs(GlobalVector& v) const = 0;

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

    virtual void setReleaseNodalForces(GlobalVector const* /*r_neq*/,
                                       int const /*process_id*/) = 0;

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
    virtual void applyKnownSolutionsPicard(
        GlobalMatrix& A, GlobalVector& rhs, GlobalVector& x,
        MathLib::DirichletBCApplicationMode const mode) const = 0;

    //! Returns whether the assembled matrix \f$A\f$ has changed and the linear
    //! solver must perform the MathLib::EigenLinearSolver::compute() step.
    virtual MathLib::LinearSolverBehaviour linearSolverNeedsToCompute()
        const = 0;
};

//! @}
}  // namespace NumLib
