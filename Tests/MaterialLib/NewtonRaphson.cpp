/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>
#include <limits>

#include "MaterialLib/SolidModels/NewtonRaphson.h"

TEST(MaterialLibNewtonRaphson, Sqrt3)
{
    static const int N = 1;  // Problem's size.

    using LocalJacobianMatrix = Eigen::Matrix<double, N, N, Eigen::RowMajor>;
    using LocalResidualVector = Eigen::Matrix<double, N, 1>;

    Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(N);
    LocalJacobianMatrix jacobian;

    const int maximum_iterations = 5;
    const double tolerance = 1e-15;

    // Solve f(x) = x^2 - S == 0, where S = 3.
    //
    // Initial value
    double state = 1;

    auto const update_jacobian = [&state](LocalJacobianMatrix& jacobian) {
        jacobian(0, 0) = 2 * state;
    };

    auto const update_residual = [&state](LocalResidualVector& residual) {
        residual[0] = state * state - 3;
    };

    auto const update_solution = [&state](
        LocalResidualVector const& increment) { state += increment[0]; };

    auto const newton_solver = MaterialLib::Solids::NewtonRaphson<
        decltype(linear_solver), LocalJacobianMatrix, decltype(update_jacobian),
        LocalResidualVector, decltype(update_residual),
        decltype(update_solution)>(linear_solver, update_jacobian,
                                   update_residual, update_solution,
                                   maximum_iterations, tolerance);
    auto const success_iterations = newton_solver.solve(jacobian);

    EXPECT_TRUE(static_cast<bool>(success_iterations));
    EXPECT_LE(*success_iterations, maximum_iterations);
    ASSERT_LE(state - std::sqrt(3), std::numeric_limits<double>::epsilon());
}
