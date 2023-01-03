/**
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "TestTools.h"

TEST(UnitTestUtilsGTestPredicateEigen, PredicateSucceeds)
{
    Tests::EigenIsNear const pred{};

    Eigen::Matrix2d const A = Eigen::Matrix2d::Random();
    EXPECT_PRED_FORMAT2(pred, A, A);

    // same test with custom error tolerance
    EXPECT_PRED_FORMAT3(pred, A, A, std::numeric_limits<double>::epsilon());
}

TEST(UnitTestUtilsGTestPredicateEigen, PredicateFailsDueToDifferentValues)
{
    Tests::EigenIsNear const pred{};

    Eigen::Matrix3d const A = Eigen::Matrix3d::Random();

    // test with custom error tolerance
    Eigen::Matrix3d const B = A + Eigen::Matrix3d::Constant(1);
    EXPECT_FALSE(pred("A", "B", "epsilon", A, B, 0.5))
        << "Predicate is wrong, because B differs more than 0.5 from A";

    // test with the default error tolerance
    Eigen::Matrix3d const C =
        A +
        Eigen::Matrix3d::Constant(2 * std::numeric_limits<double>::epsilon());
    EXPECT_FALSE(pred("A", "C", A, C))
        << "Predicate is wrong, because C differs more than the default "
           "tolerance from A";

    // only perturb a single matrix entry
    Eigen::Matrix3d D = A;
    D(2, 2) += 0.6;
    EXPECT_FALSE(pred("A", "D", "epsilon", A, D, 0.5))
        << "Predicate is wrong, because D differs more than 0.5 from A";
}

TEST(UnitTestUtilsGTestPredicateEigen, PredicateFailsDueToDifferentDimensions)
{
    Tests::EigenIsNear const pred{};

    Eigen::MatrixXd const A = Eigen::MatrixXd::Random(3, 1);
    Eigen::Vector2d const b = Eigen::Vector2d::Random();

    EXPECT_FALSE(pred("A", "b", A, b))
        << "A's and b's dimensions do not match.";

    // same test with custom error tolerance (which does not matter)
    EXPECT_FALSE(pred("A", "b", "eps", A, b, 0.5))
        << "A's and b's dimensions do not match.";
}
