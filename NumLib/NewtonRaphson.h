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

#include <Eigen/Dense>
#include <optional>

#include "BaseLib/Logging.h"

namespace NumLib
{

struct NewtonRaphsonSolverParameters
{
    int maximum_iterations;
    double residuum_tolerance;
    double increment_tolerance;
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
          _residuum_tolerance_squared(solver_parameters.residuum_tolerance *
                                      solver_parameters.residuum_tolerance),
          _increment_tolerance_squared(solver_parameters.increment_tolerance *
                                       solver_parameters.increment_tolerance)
    {
    }

    /// Returns true and the iteration number if succeeded, otherwise false and
    /// an undefined value for the number of iterations.
    std::optional<int> solve(JacobianMatrix& jacobian) const
    {
        int iteration = 0;
        ResidualVector increment;
        ResidualVector residual;
        do
        {
            // The jacobian and the residual are updated simultaneously to keep
            // consistency. The jacobian is used after the non-linear solver
            // onward.
            _jacobian_update(jacobian);
            _residual_update(residual);

            if (residual.squaredNorm() < _residuum_tolerance_squared)
            {
                break;  // convergence criteria fulfilled.
            }

            increment.noalias() =
                _linear_solver.compute(jacobian).solve(-residual);
            // DBUG("Local linear solver accuracy |J dx - r| = {:g}",
            //      (jacobian * increment + residual).norm());

            _solution_update(increment);

            if (increment.squaredNorm() < _increment_tolerance_squared)
            {
                break;  // increment to small.
            }

            // DBUG("Local Newton: Iteration #{:d} |dx| = {:g}, |r| = {:g}",
            //      iteration, increment.norm(), residual.norm());
            // fmt::print("Local Newton: Increment {}\n",
            //      fmt::join(increment.data(),
            //                increment.data() + increment.size(), ", "));
            // fmt::print("Local Newton: Residuum {}\n",
            //      fmt::join(residual.data(),
            //                residual.data() + residual.size(), ", "));
        } while (iteration++ < _maximum_iterations);

        if (iteration > _maximum_iterations)
        {
            ERR("The local Newton method did not converge within the given "
                "number of iterations. Iteration: {:d}, increment {:g}, "
                "residual: "
                "{:g}",
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
    const int _maximum_iterations;  ///< Maximum number of iterations.
    const double
        _residuum_tolerance_squared;  ///< Error tolerance for the residuum.
    const double
        _increment_tolerance_squared;  ///< Error tolerance for the increment.
};
}  // namespace NumLib
