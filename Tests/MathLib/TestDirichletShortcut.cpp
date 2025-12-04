/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <random>
#include <range/v3/numeric/iota.hpp>

#ifndef USE_PETSC
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenTools.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"

TEST(MathLib, DirichletShortcut)
{
    Eigen::Index const nrows = 100;
    Eigen::Index const ncols = 100;
    Eigen::Index const sparsity = 10;

    std::vector<Eigen::Index> ids(ncols);
    ranges::iota(ids, 0);

    std::mt19937 gen{std::random_device{}()};
    std::uniform_real_distribution dist{1., 2.};

    // fill A
    std::vector<Eigen::Index> selected_col_ids;

    MathLib::EigenMatrix A{nrows, sparsity};

    for (Eigen::Index row = 0; row < nrows; ++row)
    {
        selected_col_ids.clear();
        std::ranges::sample(ids, std::back_inserter(selected_col_ids), sparsity,
                            gen);

        for (auto col : selected_col_ids)
        {
            A.setValue(row, col, dist(gen));
        }
    }

    // fill rhs
    MathLib::EigenVector rhs{nrows};
    for (Eigen::Index row = 0; row < nrows; ++row)
    {
        rhs[row] = dist(gen);
    }

    // fill ids and known values
    std::vector<Eigen::Index> known_x_ids;
    std::ranges::sample(ids, std::back_inserter(known_x_ids), sparsity, gen);
    std::vector<double> known_x_values(known_x_ids.size());
    for (auto& known_x : known_x_values)
    {
        known_x = dist(gen);
    }

    auto A_copy1 = A;
    auto rhs_copy1 = rhs;
    MathLib::EigenVector x{nrows};

    // reference
    MathLib::applyKnownSolution(
        A_copy1, rhs_copy1, x, known_x_ids, known_x_values,
        MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE);

    auto A_copy2 = A;
    auto rhs_copy2 = rhs;

    // testee
    MathLib::applyKnownSolution(
        A_copy2, rhs_copy2, x, known_x_ids, known_x_values,
        MathLib::DirichletBCApplicationMode::FAST_INCOMPLETE_MATRIX_UPDATE);

    EXPECT_THAT(
        rhs_copy2.getRawVector(),
        testing::Pointwise(
            testing::DoubleNear(20 * std::numeric_limits<double>::epsilon()),
            rhs_copy1.getRawVector()));
}
#endif
