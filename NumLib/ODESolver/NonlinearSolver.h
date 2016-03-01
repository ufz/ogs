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

#include "MathLib/LinAlg/LinearSolver.h"

#include "Types.h"
#include "NonlinearSystem.h"

namespace BaseLib
{
class ConfigTree;
}

// TODO Document in the ODE solver lib, which matrices and vectors that are passed around
//      as method arguments are guaranteed to be of the right size (and zeroed out) and
//      which are not.

namespace NumLib
{

/*! Common interface for nonlinear solvers.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template<typename Matrix, typename Vector>
class NonlinearSolverBase
{
public:
    /*! Only assemble the equation system.
     *
     * \note This method is needed to preload CrankNicolson time discretization scheme.
     *       It is not used for the general solver steps; in those only the solve() method
     *       is needed.
     *
     * \param x   the state at which the equation system will be assembled.
     */
    virtual void assemble(Vector const& x) const = 0;

    /*! Assemble and solve the equation system.
     *
     * \param x   in: the initial guess, out: the solution.
     *
     * \retval true if the equation system could be solved
     * \retval false otherwise
     */
    virtual bool solve(Vector& x) = 0;

    virtual ~NonlinearSolverBase() = default;
};

//! \addtogroup ODESolver
//! @{

/*! Find a solution to a nonlinear equation.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 * \tparam NLTag  a tag indicating the method used for solving the equation.
 */
template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class NonlinearSolver;


/*! Find a solution to a nonlinear equation using the Newton-Raphson method.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template<typename Matrix, typename Vector>
class NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Newton> final
        : public NonlinearSolverBase<Matrix, Vector>
{
public:
    //! Type of the nonlinear equation system to be solved.
    using System = NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton>;
    using LinearSolver = MathLib::LinearSolver<Matrix, Vector>;

    /*! Constructs a new instance.
     *
     * \param linear_solver the linear solver used by this nonlinear solver.
     * \param tol     the tolerance of the solver. \todo Be more specific about that!
     * \param maxiter the maximum number of iterations used to solve the equation.
     */
    explicit
    NonlinearSolver(LinearSolver& linear_solver,
                    double const tol, const unsigned maxiter)
        : _linear_solver(linear_solver)
        , _tol(tol)
        , _maxiter(maxiter)
    {}

    //! Set the nonlinear equation system that will be solved.
    void setEquationSystem(System& eq) { _equation_system = &eq; }

    void assemble(Vector const& x) const override;

    bool solve(Vector& x) override;

private:
    LinearSolver& _linear_solver;
    System*       _equation_system = nullptr;

    const double _tol;       //!< tolerance of the solver
    const unsigned _maxiter; //!< maximum number of iterations

    Vector _res;             //!< The residual.
    Matrix _J;               //!< The Jacobian of the residual.
    Vector _minus_delta_x;   //!< The Newton-Raphson method solves the linearized equation
                             //!< \f$ J \cdot (-\Delta x) = r \f$ repeatedly.
    double const _alpha = 1; //!< Damping factor. \todo Add constructor parameter.
};


/*! Find a solution to a nonlinear equation using the Picard fixpoint iteration method.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the equation.
 */
template<typename Matrix, typename Vector>
class NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Picard> final
        : public NonlinearSolverBase<Matrix, Vector>
{
public:
    //! Type of the nonlinear equation system to be solved.
    using System = NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard>;
    using LinearSolver = MathLib::LinearSolver<Matrix, Vector>;

    /*! Constructs a new instance.
     *
     * \param linear_solver the linear solver used by this nonlinear solver.
     * \param tol     the tolerance of the solver. \todo Be more specific about that!
     * \param maxiter the maximum number of iterations used to solve the equation.
     */
    explicit
    NonlinearSolver(LinearSolver& linear_solver,
                    double const tol, const unsigned maxiter)
        : _linear_solver(linear_solver)
        , _tol(tol)
        , _maxiter(maxiter)
    {}

    //! Set the nonlinear equation system that will be solved.
    void setEquationSystem(System& eq) { _equation_system = &eq; }

    void assemble(Vector const& x) const override;

    bool solve(Vector& x) override;

private:
    LinearSolver& _linear_solver;
    System*       _equation_system = nullptr;

    const double _tol;       //!< tolerance of the solver
    const unsigned _maxiter; //!< maximum number of iterations

    Matrix _A;     //!< \c Matrix describing the linearized system.
    Vector _rhs;   //!< \c Vector describing the linearized system.
    Vector _x_new; //!< \c Vector to store solutions of \f$ A \cdot x = \mathit{rhs} \f$.
};

/*! Creates a new nonlinear solver from the given configuration.
 *
 * \param linear_solver the linear solver that will be used by the nonlinear
 *        solver
 * \param config configuration settings
 *
 * \return a pair <tt>(nl_slv, tag)</tt> where \c nl_slv is the generated
 *         nonlinear solver instance and the \c tag indicates if it uses the
 *         Picard or Newton-Raphson method
 */
template<typename Matrix, typename Vector>
std::pair<
    std::unique_ptr<NonlinearSolverBase<Matrix, Vector> >,
    NonlinearSolverTag
>
createNonlinearSolver(MathLib::LinearSolver<Matrix, Vector>& linear_solver,
                      BaseLib::ConfigTree const& config);

//! @}

} // namespace NumLib

#include "NonlinearSolver-impl.h"

#endif // NUMLIB_NONLINEARSOLVER_H
