/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 23, 2024
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <random>

#include "MathLib/Nonlinear/CubicRoots.h"

// Test a simple cubic equation: x^3 - 6x^2 + 11x - 6 = 0
// This equation has known roots: x = 1, 2, 3
TEST(MathLibCubicSolver, SimpleCubicEquation)
{
    MathLib::CubicSolver solver(1.0, -6.0, 11.0, -6.0);
    auto roots = solver.solve();

    ASSERT_EQ(roots.size(), 3);  // There should be 3 real roots

    // Check if the known roots are in the result
    EXPECT_NEAR(roots[0], 1.0, 1e-6);
    EXPECT_NEAR(roots[1], 2.0, 1e-6);
    EXPECT_NEAR(roots[2], 3.0, 1e-6);
}

// Test a cubic equation with one real root: x^3 - 3x^2 + 3x - 1 = 0
// This equation has a triple root at x = 1
TEST(MathLibCubicSolver, CubicEquationWithTripleRoot)
{
    MathLib::CubicSolver solver(1.0, -3.0, 3.0, -1.0);
    auto roots = solver.solve();

    ASSERT_EQ(roots.size(), 3);  // There should be 3 roots in the result

    // Check if all three roots are approximately 1.0
    EXPECT_NEAR(roots[0], 1.0, 1e-6);
    EXPECT_NEAR(roots[1], 1.0, 1e-6);
    EXPECT_NEAR(roots[2], 1.0, 1e-6);
}

// Test a cubic equation with one real root: x^3 + 7x^2 - 12x - 30 = 0
// This equation has two negative and one positive reel at x = 2.515
TEST(MathLibCubicSolver, SmallestPositiveRoot)
{
    MathLib::CubicSolver solver(1.0, 7.0, -12.0, -30.0);
    auto roots = solver.solve();
    auto positive_root = solver.smallestPositiveRealRoot();

    ASSERT_EQ(roots.size(), 3);  // There should be 3 roots in the result

    // Check if the smallest positive root is about 2.515
    EXPECT_NEAR(positive_root, 2.5148929484550, 1e-6);
}

// Test a cubic equation with one real root and complex roots
// x^3 - x + 1 = 0 has one real root and two complex roots
TEST(MathLibCubicSolver, CubicEquationWithComplexRoots)
{
    MathLib::CubicSolver solver(1.0, 0.0, -1.0, 1.0);
    auto roots = solver.solve();

    ASSERT_EQ(roots.size(), 3);  // Only one real root should be found
    EXPECT_NEAR(roots[0], -1.3247,
                1e-4);  // Check against the expected real root
}

// Test an equation where 'a' is effectively zero, which should throw an
// exception
TEST(MathLibCubicSolver, DegenerateCase)
{
    EXPECT_THROW(MathLib::CubicSolver solver(0.0, -6.0, 11.0, -6.0),
                 std::runtime_error);
}

// Helper function to evaluate the polynomial ax^3 + bx^2 + cx + d at a given x
double evaluatePolynomial(double a, double b, double c, double d, double x)
{
    return a * x * x * x + b * x * x + c * x + d;
}

// Helper function to perform the root finding test on a given set of roots
void performCubicSolverTest(
    const std::vector<std::tuple<double, double, double>>& root_sets,
    const double tolerance = 1e-3)
{
    for (const auto& roots_tuple : root_sets)
    {
        double r1, r2, r3;
        std::tie(r1, r2, r3) = roots_tuple;

        // Construct the coefficients of the polynomial using the given formula
        const double a = 1.0;  // Coefficient of x^3
        const double b = -(r1 + r2 + r3);
        const double c = (r1 * r2 + r1 * r3 + r2 * r3);
        const double d = -r1 * r2 * r3;

        // Instantiate the solver
        MathLib::CubicSolver solver(a, b, c, d);
        auto roots = solver.solve();

        // There should always be 3 roots (even if some are repeated)
        ASSERT_EQ(roots.size(), 3);

        // Check each root by substituting it back into the polynomial
        for (const auto& root : roots)
        {
            const double value = evaluatePolynomial(a, b, c, d, root);

            // Check if the polynomial evaluated at the root is approximately
            // zero
            EXPECT_NEAR(value, 0.0, tolerance)
                << "Root " << root
                << " did not satisfy the equation for polynomial (" << a
                << "x^3 + " << b << "x^2 + " << c << "x + " << d << ")";
        }
    }
}

// Test with predefined root sets
TEST(MathLibCubicSolver, DefinedCubicEquations)
{
    // Define a vector of tuples representing sets of desired roots (r1, r2, r3)
    std::vector<std::tuple<double, double, double>> defined_root_sets = {
        {1.0, 2.0, 3.0},     // Example 1: Three distinct roots
        {1.0, 1.0, 2.0},     // Example 2: One double root and one single root
        {1.0, 1.0, 1.0},     // Example 3: Triple root
        {0.0, 2.0, -2.0},    // Example 4: Includes zero root
        {-1.0, -1.0, 3.0},   // Example 5: Includes negative roots
        {-1.0, -5.0, -3.0},  // Example 6: Unsorted roots
        {-100.0, -5000.0, -300.0},  // Example 7: Large range of roots
        {100.0, 500.0, 3000.0},
    };

    // Use the helper function to perform the test
    performCubicSolverTest(defined_root_sets);
}

// Test with randomly generated root sets
TEST(MathLibCubicSolver, RandomCubicEquations)
{
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(
        0.0,
        2000.0);  // The roots of the polynomial are within this distribution.

    std::vector<std::tuple<double, double, double>> random_root_sets;

    // Generate 100k sets of random roots
    for (int i = 0; i < 100000; ++i)
    {
        double r1 = dis(gen);
        double r2 = dis(gen);
        double r3 = dis(gen);

        // Sort the roots in ascending order
        std::array<double, 3> roots = {r1, r2, r3};
        std::sort(roots.begin(), roots.end());

        random_root_sets.emplace_back(roots[0], roots[1], roots[2]);
    }

    // Use the helper function to perform the test
    performCubicSolverTest(random_root_sets);
}