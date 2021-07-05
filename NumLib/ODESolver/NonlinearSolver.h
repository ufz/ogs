/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include "ConvergenceCriterion.h"
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
        std::function<void(int, std::vector<GlobalVector*> const&)> const&
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
     * \param linear_solver the linear solver used by this nonlinear solver.
     * \param maxiter the maximum number of iterations used to solve the
     *                equation.
     * \param damping A positive damping factor.
     * \see _damping
     */
    explicit NonlinearSolver(GlobalLinearSolver& linear_solver,
                             int const maxiter,
                             double const damping = 1.0)
        : _linear_solver(linear_solver), _maxiter(maxiter), _damping(damping)
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
        std::function<void(int, std::vector<GlobalVector*> const&)> const&
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
    int const _maxiter;  //!< maximum number of iterations

    //! A positive damping factor. The default value 1.0 gives a non-damped
    //! Newton method. Common values are in the range 0.5 to 0.7 for somewhat
    //! conservative method and seldom become smaller than 0.2 for very
    //! conservative approach.
    double const _damping;

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
        std::function<void(int, std::vector<GlobalVector*> const&)> const&
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

/*! Creates a new nonlinear solver from the given configuration.
 *
 * \param linear_solver the linear solver that will be used by the nonlinear
 *                      solver
 * \param config configuration settings
 *
 * \return a pair <tt>(nl_slv, tag)</tt> where \c nl_slv is the generated
 *         nonlinear solver instance and the \c tag indicates if it uses
 *         the Picard or Newton-Raphson method
 */
std::pair<std::unique_ptr<NonlinearSolverBase>, NonlinearSolverTag>
createNonlinearSolver(GlobalLinearSolver& linear_solver,
                      BaseLib::ConfigTree const& config);

//! @}

}  // namespace NumLib
