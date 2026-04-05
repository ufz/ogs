// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <utility>

#include "ConvergenceCriterion.h"
#include "MathLib/LinAlg/GlobalLinearSolverType.h"
#include "NewtonStepStrategy.h"
#include "NonlinearSolverStatus.h"
#include "NonlinearSystem.h"
#include "Types.h"

namespace BaseLib
{
class ConfigTree;
}

// TODO Document in the ODE solver lib, which matrices and vectors that are
// passed around as method arguments are guaranteed to be of the right size
// (and zeroed out) and which are not.

namespace NumLib
{
/*! Common interface for nonlinear solvers.
 *
 */
class NonlinearSolverBase
{
public:
    virtual void calculateNonEquilibriumInitialResiduum(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id) = 0;

    /*! Assemble and solve the equation system.
     *
     * \param x   in: the initial guess, out: the solution.
     * \param x_prev previous time step solution.
     * \param postIterationCallback called after each iteration if set.
     * \param process_id usually used in staggered schemes.
     *
     * \retval true if the equation system could be solved
     * \retval false otherwise
     */
    virtual NonlinearSolverStatus solve(
        std::vector<GlobalVector*>& x,
        std::vector<GlobalVector*> const& x_prev,
        std::function<void(int, bool, std::vector<GlobalVector*> const&)> const&
            postIterationCallback,
        int const process_id) = 0;

    virtual ~NonlinearSolverBase() = default;
};

//! \addtogroup ODESolver
//! @{

/*! Find a solution to a nonlinear equation.
 *
 * \tparam NLTag  a tag indicating the method used for solving the equation.
 */
template <NonlinearSolverTag NLTag>
class NonlinearSolver;

/*! Find a solution to a nonlinear equation using the Newton-Raphson method.
 */
template <>
class NonlinearSolver<NonlinearSolverTag::Newton> final
    : public NonlinearSolverBase
{
public:
    //! Type of the nonlinear equation system to be solved.
    using System = NonlinearSystem<NonlinearSolverTag::Newton>;

    /*! Constructs a new instance.
     *
     * \param linear_solver    the linear solver used by this nonlinear solver.
     * \param maxiter          the maximum number of iterations used to solve
     *                         the equation.
     * \param newton_strategy  globalization / step-acceptance strategy
     *                         (e.g. fixed damping, line search).  Ownership
     *                         is transferred to the solver.
     * \param recompute_jacobian recompute the Jacobian every this many steps.
     */
    explicit NonlinearSolver(GlobalLinearSolver& linear_solver,
                             int const maxiter,
                             std::unique_ptr<NewtonStepStrategy>
                                 newton_strategy,
                             int const recompute_jacobian = 1)
        : _linear_solver(linear_solver),
          _maxiter(maxiter),
          _step_strategy(std::move(newton_strategy)),
          _recompute_jacobian(recompute_jacobian)
    {
    }

    ~NonlinearSolver();

    //! Set the nonlinear equation system that will be solved and the
    //! convergence criterion used to terminate the iteration.  Also forwards
    //! the criterion to the step strategy so that strategies that depend on
    //! residual or delta-x information (e.g. non-negative damping) can use it.
    void setEquationSystem(System& eq, ConvergenceCriterion& conv_crit)
    {
        _equation_system = &eq;
        _convergence_criterion = &conv_crit;
        _step_strategy->setConvergenceCriterion(*_convergence_criterion);
    }

    void calculateNonEquilibriumInitialResiduum(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev,
        int const process_id) override;

    NonlinearSolverStatus solve(
        std::vector<GlobalVector*>& x,
        std::vector<GlobalVector*> const& x_prev,
        std::function<void(int, bool, std::vector<GlobalVector*> const&)> const&
            postIterationCallback,
        int const process_id) override;

    void compensateNonEquilibriumInitialResiduum(bool const value)
    {
        _compensate_non_equilibrium_initial_residuum = value;
    }

private:
    GlobalLinearSolver& _linear_solver;
    System* _equation_system = nullptr;

    int const _maxiter;  //!< maximum number of iterations

    //! Globalization / step-acceptance strategy (e.g. fixed damping).
    std::unique_ptr<NewtonStepStrategy> _step_strategy;

    //! Convergence criterion used to terminate the Newton iteration.
    ConvergenceCriterion* _convergence_criterion = nullptr;

    int const _recompute_jacobian =
        1;  //!< Recompute Jacobian every this many steps.

    GlobalVector* _r_neq = nullptr;      //!< non-equilibrium initial residuum.
    std::size_t _res_id = 0u;            //!< ID of the residual vector.
    std::size_t _J_id = 0u;              //!< ID of the Jacobian matrix.
    std::size_t _minus_delta_x_id = 0u;  //!< ID of the \f$ -\Delta x\f$ vector.
    std::size_t _x_new_id =
        0u;  //!< ID of the vector storing \f$ x - (-\Delta x) \f$.
    std::size_t _r_neq_id = 0u;  //!< ID of the non-equilibrium initial
                                 //! residuum vector.

    /// Enables computation of the non-equilibrium initial residuum \f$ r_{\rm
    /// neq} \f$ before the first time step. The forces are zero if the external
    /// forces are in equilibrium with the initial state/initial conditions.
    /// During the simulation the new residuum reads \f$ \tilde r = r - r_{\rm
    /// neq} \f$.
    bool _compensate_non_equilibrium_initial_residuum = false;
};

/*! Find a solution to a nonlinear equation using the Picard fixpoint iteration
 * method.
 *
 */
template <>
class NonlinearSolver<NonlinearSolverTag::Picard> final
    : public NonlinearSolverBase
{
public:
    //! Type of the nonlinear equation system to be solved.
    using System = NonlinearSystem<NonlinearSolverTag::Picard>;

    /*! Constructs a new instance.
     *
     * \param linear_solver the linear solver used by this nonlinear solver.
     * \param maxiter the maximum number of iterations used to solve the
     *                equation.
     */
    explicit NonlinearSolver(GlobalLinearSolver& linear_solver,
                             const int maxiter)
        : _linear_solver(linear_solver), _maxiter(maxiter)
    {
    }

    ~NonlinearSolver();

    //! Set the nonlinear equation system that will be solved.
    //! TODO doc
    void setEquationSystem(System& eq, ConvergenceCriterion& conv_crit)
    {
        _equation_system = &eq;
        _convergence_criterion = &conv_crit;
    }

    void calculateNonEquilibriumInitialResiduum(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev,
        int const process_id) override;

    NonlinearSolverStatus solve(
        std::vector<GlobalVector*>& x,
        std::vector<GlobalVector*> const& x_prev,
        std::function<void(int, bool, std::vector<GlobalVector*> const&)> const&
            postIterationCallback,
        int const process_id) override;

    void compensateNonEquilibriumInitialResiduum(bool const value)
    {
        _compensate_non_equilibrium_initial_residuum = value;
    }

private:
    GlobalLinearSolver& _linear_solver;
    System* _equation_system = nullptr;

    // TODO doc
    ConvergenceCriterion* _convergence_criterion = nullptr;
    const int _maxiter;  //!< maximum number of iterations

    GlobalVector* _r_neq = nullptr;  //!< non-equilibrium initial residuum.
    std::size_t _A_id = 0u;          //!< ID of the \f$ A \f$ matrix.
    std::size_t _rhs_id = 0u;        //!< ID of the right-hand side vector.
    std::size_t _x_new_id = 0u;  //!< ID of the vector storing the solution of
                                 //! the linearized equation.
    std::size_t _r_neq_id = 0u;  //!< ID of the non-equilibrium initial
                                 //! residuum vector.

    // clang-format off
    /// \copydoc NumLib::NonlinearSolver<NonlinearSolverTag::Newton>::_compensate_non_equilibrium_initial_residuum
    bool _compensate_non_equilibrium_initial_residuum = false;
    // clang-format on
};
//! @}

}  // namespace NumLib
