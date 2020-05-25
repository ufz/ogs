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

#include <boost/optional.hpp>
#include "BaseLib/Logging.h"

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
        : linear_solver_(linear_solver),
          jacobian_update_(jacobian_update),
          residual_update_(residual_update),
          solution_update_(solution_update),
          maximum_iterations_(solver_parameters.maximum_iterations),
          tolerance_squared_(solver_parameters.error_tolerance *
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
            jacobian_update_(jacobian);
            residual_update_(residual);

            if (residual.squaredNorm() < tolerance_squared_)
            {
                break;  // convergence criteria fulfilled.
            }

            increment.noalias() =
                linear_solver_.compute(jacobian).solve(-residual);
            // DBUG("Local linear solver accuracy |J dx - r| = {:g}",
            //      (jacobian * increment + residual).norm());

            solution_update_(increment);

            // DBUG("Local Newton: Iteration #{:d} |dx| = {:g}, |r| = {:g}",
            //      iteration, increment.norm(), residual.norm());
        } while (iteration++ < maximum_iterations_);

        if (iteration > maximum_iterations_)
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
    LinearSolver& linear_solver_;
    JacobianMatrixUpdate jacobian_update_;
    ResidualUpdate residual_update_;
    SolutionUpdate solution_update_;
    const int maximum_iterations_;  ///< Maximum number of iterations.
    const double tolerance_squared_;    ///< Error tolerance for the residual.
};
}  // namespace NumLib
