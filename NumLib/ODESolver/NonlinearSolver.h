/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_NONLINEARSOLVER_H
#define NUMLIB_NONLINEARSOLVER_H

#include <memory>
#include <utility>
#include <logog/include/logog.hpp>

#include "ConvergenceCriterion.h"
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
    /*! Only assemble the equation system.
     *
     * \note This method is needed to preload CrankNicolson time discretization
     * scheme. It is not used for the general solver steps; in those only the
     * solve() method is needed.
     *
     * \param x   the state at which the equation system will be assembled.
     */
    virtual void assemble(GlobalVector const& x) const = 0;

    /*! Assemble and solve the equation system.
     *
     * \param x   in: the initial guess, out: the solution.
     * \param postIterationCallback called after each iteration if set.
     *
     * \retval true if the equation system could be solved
     * \retval false otherwise
     */
    virtual bool solve(GlobalVector& x,
                       std::function<void(unsigned, GlobalVector const&)> const&
                           postIterationCallback) = 0;

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
 *
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
     */
    explicit NonlinearSolver(
        GlobalLinearSolver& linear_solver,
        const unsigned maxiter)
        : _linear_solver(linear_solver),
          _maxiter(maxiter)
    {
    }

    //! Set the nonlinear equation system that will be solved.
    //! TODO doc
    void setEquationSystem(System& eq, ConvergenceCriterion& conv_crit)
    {
        _equation_system = &eq;
        _convergence_criterion = &conv_crit;
    }
    void assemble(GlobalVector const& x) const override;

    bool solve(GlobalVector& x,
               std::function<void(unsigned, GlobalVector const&)> const&
                   postIterationCallback) override;

private:
    GlobalLinearSolver& _linear_solver;
    System* _equation_system = nullptr;

    // TODO doc
    ConvergenceCriterion* _convergence_criterion = nullptr;
    const unsigned _maxiter;  //!< maximum number of iterations

    double const _alpha =
        1;  //!< Damping factor. \todo Add constructor parameter.

    std::size_t _res_id = 0u;            //!< ID of the residual vector.
    std::size_t _J_id = 0u;              //!< ID of the Jacobian matrix.
    std::size_t _minus_delta_x_id = 0u;  //!< ID of the \f$ -\Delta x\f$ vector.
    std::size_t _x_new_id =
        0u;  //!< ID of the vector storing \f$ x - (-\Delta x) \f$.
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
                             const unsigned maxiter)
        : _linear_solver(linear_solver), _maxiter(maxiter)
    {
    }

    //! Set the nonlinear equation system that will be solved.
    //! TODO doc
    void setEquationSystem(System& eq, ConvergenceCriterion& conv_crit)
    {
        _equation_system = &eq;
        _convergence_criterion = &conv_crit;
    }
    void assemble(GlobalVector const& x) const override;

    bool solve(GlobalVector& x,
               std::function<void(unsigned, GlobalVector const&)> const&
                   postIterationCallback) override;

private:
    GlobalLinearSolver& _linear_solver;
    System* _equation_system = nullptr;

    // TODO doc
    ConvergenceCriterion* _convergence_criterion = nullptr;
    const unsigned _maxiter;  //!< maximum number of iterations

    std::size_t _A_id = 0u;      //!< ID of the \f$ A \f$ matrix.
    std::size_t _rhs_id = 0u;    //!< ID of the right-hand side vector.
    std::size_t _x_new_id = 0u;  //!< ID of the vector storing the solution of
                                 //!the linearized equation.
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
createNonlinearSolver(
    GlobalLinearSolver& linear_solver,
    BaseLib::ConfigTree const& config);

//! @}

}  // namespace NumLib

#endif  // NUMLIB_NONLINEARSOLVER_H
