/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/optional.hpp>
#include <logog/include/logog.hpp>

#include <Eigen/Dense>

namespace NumLib
{

struct NewtonRaphsonSolverParameters
{
    int const maximum_iterations;
    double const error_tolerance;
};

/// Newton-Raphson solver for system of equations using an Eigen linear solvers
/// library.
/// The current implementation does not update the solution itself, but calls a
/// function for the solution's update with the current increment.
template <typename LinearSolver, typename JacobianMatrix,
          typename JacobianMatrixUpdate, typename ResidualVector,
          typename ResidualUpdate, typename SolutionUpdate>
class NewtonRaphson final
{
public:
    NewtonRaphson(LinearSolver& linear_solver,
                  JacobianMatrixUpdate jacobian_update,
                  ResidualUpdate residual_update,
                  SolutionUpdate solution_update,
                  NewtonRaphsonSolverParameters const& solver_parameters)
        : _linear_solver(linear_solver),
          _jacobian_update(jacobian_update),
          _residual_update(residual_update),
          _solution_update(solution_update),
          _maximum_iterations(solver_parameters.maximum_iterations),
          _tolerance_squared(solver_parameters.error_tolerance *
                             solver_parameters.error_tolerance)
    {
    }

    /// Returns true and the iteration number if succeeded, otherwise false and
    /// an undefined value for the number of iterations.
    boost::optional<int> solve(JacobianMatrix& jacobian) const
    {
        int iteration = 0;
        ResidualVector increment;
        ResidualVector residual;
        do
        {
            // The jacobian and the residual are updated simulataniously to keep
            // consistency. The jacobian is used after the non-linear solver
            // onward.
            _jacobian_update(jacobian);
            _residual_update(residual);

            if (residual.squaredNorm() < _tolerance_squared)
                break;  // convergence criteria fulfilled.

            increment.noalias() =
                _linear_solver.compute(jacobian).solve(-residual);
            // DBUG("Local linear solver accuracy |J dx - r| = %g",
            //      (jacobian * increment + residual).norm());

            _solution_update(increment);

            // DBUG("Local Newton: Iteration #%d |dx| = %g, |r| = %g",
            //      iteration, increment.norm(), residual.norm());
        } while (iteration++ < _maximum_iterations);

        if (iteration > _maximum_iterations)
        {
            ERR("The local Newton method did not converge within the given "
                "number of iterations. Iteration: %d, increment %g, residual: "
                "%g",
                iteration - 1, increment.norm(), residual.norm());
            return {};
        }

        return iteration;
    };

private:
    LinearSolver& _linear_solver;
    JacobianMatrixUpdate _jacobian_update;
    ResidualUpdate _residual_update;
    SolutionUpdate _solution_update;
    const int _maximum_iterations;
    const double _tolerance_squared;
};
}  // namespace NumLib
