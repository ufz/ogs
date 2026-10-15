// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <random>
#include <range/v3/algorithm/sample.hpp>
#include <range/v3/numeric/iota.hpp>
#include <string>
#include <vector>

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

namespace
{
/// The former implementation of the complete matrix update, which transposes
/// the matrix twice, as reference. Its right-hand-side update is written as
/// an explicit std::fma, the rounding the former "b -= a*x" gets whenever the
/// compiler contracts it, so the reference's rounding does not depend on the
/// contraction policy of the build. The new implementation must reproduce it
/// exactly, since benchmarks with tight tolerances depend on it.
void applyKnownSolutionCompleteByTransposition(
    MathLib::EigenMatrix::RawMatrixType& A,
    MathLib::EigenVector::RawVectorType& b,
    std::vector<MathLib::EigenMatrix::IndexType> const& ids,
    std::vector<double> const& values)
{
    using SpMat = MathLib::EigenMatrix::RawMatrixType;

    if (ids.size() != values.size())
    {
        OGS_FATAL("Reference: ids and values must pair up: {} ids, {} values.",
                  ids.size(), values.size());
    }

    for (auto const row_id : ids)
    {
        for (SpMat::InnerIterator it(A, row_id); it; ++it)
        {
            if (it.col() != row_id)
            {
                it.valueRef() = 0.0;
            }
        }
    }

    SpMat AT = A.transpose();
    AT.reserve(Eigen::VectorXi::Constant(A.rows(), 1));

    for (std::size_t ix = 0; ix < ids.size(); ix++)
    {
        SpMat::Index const row_id = ids[ix];
        auto const x = values[ix];
        for (SpMat::InnerIterator it(AT, row_id); it; ++it)
        {
            if (it.col() == row_id)
            {
                continue;
            }
            b[it.col()] = std::fma(-it.value(), x, b[it.col()]);
            it.valueRef() = 0.0;
        }

        auto& c = AT.coeffRef(row_id, row_id);
        if (c != 0.0)
        {
            b[row_id] = x * c;
        }
        else
        {
            b[row_id] = x;
            c = 1.0;
        }
    }

    A = AT.transpose();
}
}  // namespace

TEST(MathLib, DirichletCompleteRejectsInvalidId)
{
    MathLib::EigenMatrix A{4, 3};
    A.setValue(0, 0, 1.0);
    MathLib::EigenVector b{4};
    MathLib::EigenVector x{4};

    for (auto const id :
         {MathLib::EigenMatrix::IndexType{-1}, A.getNumberOfRows()})
    {
        std::vector<MathLib::EigenMatrix::IndexType> const known_x_ids{id};
        std::vector<double> const known_x_values{1.0};
        EXPECT_THROW(
            MathLib::applyKnownSolution(
                A, b, x, known_x_ids, known_x_values,
                MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE),
            std::runtime_error);
    }
}

TEST(MathLib, DirichletCompleteWithoutKnownIdsChangesNothing)
{
    MathLib::EigenMatrix A{3};
    A.setValue(0, 0, 1.0);
    A.setValue(0, 1, 2.0);
    MathLib::EigenVector b{3};
    b[0] = 4.0;
    MathLib::EigenVector x{3};

    Eigen::MatrixXd const A_before{A.getRawMatrix()};
    MathLib::applyKnownSolution(
        A, b, x, {}, {},
        MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE);

    EXPECT_TRUE(A_before == Eigen::MatrixXd(A.getRawMatrix()));
    EXPECT_EQ(4.0, b.getRawVector()[0]);
}

TEST(MathLib, DirichletRejectsMismatchedIdsAndValues)
{
    MathLib::EigenMatrix A{4, 3};
    A.setValue(0, 0, 1.0);
    MathLib::EigenVector b{4};
    MathLib::EigenVector x{4};

    std::vector<MathLib::EigenMatrix::IndexType> const ids{0, 1};
    std::vector<MathLib::EigenMatrix::IndexType> const single_id{0};
    std::vector<double> const single_value{1.0};
    std::vector<double> const values{1.0, 2.0};

    for (auto const mode :
         {MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE,
          MathLib::DirichletBCApplicationMode::FAST_INCOMPLETE_MATRIX_UPDATE})
    {
        EXPECT_THROW(
            MathLib::applyKnownSolution(A, b, x, ids, single_value, mode),
            std::runtime_error);
        EXPECT_THROW(
            MathLib::applyKnownSolution(A, b, x, single_id, values, mode),
            std::runtime_error);
    }
}

TEST(MathLib, DirichletCompleteMatchesTransposition)
{
    Eigen::Index const n = 100;
    Eigen::Index const sparsity = 10;

    std::vector<Eigen::Index> ids(n);
    ranges::iota(ids, 0);

    std::mt19937 gen{std::random_device{}()};
    std::uniform_real_distribution dist{1., 2.};

    MathLib::EigenMatrix A{n, sparsity};
    // Row 0 stays empty as for deactivated subdomains, s.t. its diagonal
    // entry has to be inserted.
    std::vector<Eigen::Index> selected_col_ids;
    for (Eigen::Index row = 1; row < n; ++row)
    {
        selected_col_ids.clear();
        ranges::sample(ids, std::back_inserter(selected_col_ids), sparsity,
                       gen);
        for (auto col : selected_col_ids)
        {
            A.setValue(row, col, dist(gen));
        }
    }
    // An explicit zero on the diagonal of a known row.
    A.setValue(1, 1, 0.0);

    MathLib::EigenVector rhs{n};
    for (Eigen::Index row = 0; row < n; ++row)
    {
        rhs[row] = dist(gen);
    }

    // Rows 0 and 1 are known, and further distinct random ones, in random
    // order.
    constexpr std::size_t num_fixed_ids = 2;
    std::vector<Eigen::Index> const other_ids(ids.begin() + num_fixed_ids,
                                              ids.end());
    std::vector<Eigen::Index> known_x_ids{0, 1};
    ranges::sample(other_ids, std::back_inserter(known_x_ids), 2 * sparsity,
                   gen);
    std::shuffle(known_x_ids.begin(), known_x_ids.end(), gen);
    // Some ids occur twice, with different values, as for a node with a
    // Dirichlet BC that also lies on the boundary of a deactivated subdomain.
    constexpr std::size_t num_duplicated_ids = 3;
    for (std::size_t i = num_fixed_ids; i < num_fixed_ids + num_duplicated_ids;
         ++i)
    {
        known_x_ids.push_back(known_x_ids[i]);
    }
    std::vector<double> known_x_values(known_x_ids.size());
    for (auto& known_x : known_x_values)
    {
        known_x = dist(gen);
    }

    auto A_reference = A;
    auto rhs_reference = rhs;
    applyKnownSolutionCompleteByTransposition(A_reference.getRawMatrix(),
                                              rhs_reference.getRawVector(),
                                              known_x_ids, known_x_values);

    MathLib::EigenVector x{n};
    MathLib::applyKnownSolution(
        A, rhs, x, known_x_ids, known_x_values,
        MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE);

    // Same values and same sparsity pattern, including explicit zeros.
    EXPECT_TRUE(A.getRawMatrix().isCompressed());
    EXPECT_EQ(A_reference.getRawMatrix().nonZeros(),
              A.getRawMatrix().nonZeros());
    EXPECT_TRUE((Eigen::MatrixXd(A_reference.getRawMatrix()) ==
                 Eigen::MatrixXd(A.getRawMatrix())));

    // Identical right-hand side, also in round-off and for the duplicate ids.
    EXPECT_THAT(
        rhs.getRawVector(),
        testing::Pointwise(testing::Eq(), rhs_reference.getRawVector()));
    EXPECT_EQ(1.0, A.getRawMatrix().coeff(0, 0));
    EXPECT_EQ(1.0, A.getRawMatrix().coeff(1, 1));
}
/// Timing of the complete matrix update against the former implementation on
/// a matrix with the sparsity of bilinear quadrilaterals on an n x n node grid
/// (9-point stencil) and Dirichlet BCs on all boundary nodes. Disabled by
/// default; run with
///   testrunner --gtest_also_run_disabled_tests
///       --gtest_filter=MathLib.DISABLED_DirichletCompletePerformance
/// The grid size can be set via the environment variable
/// OGS_DIRICHLET_PERFORMANCE_GRID (default: 1000 nodes per direction).
TEST(MathLib, DISABLED_DirichletCompletePerformance)
{
    Eigen::Index n = 1000;
    if (char const* const grid = std::getenv("OGS_DIRICHLET_PERFORMANCE_GRID"))
    {
        try
        {
            n = std::stol(grid);
        }
        catch (std::exception const&)
        {
            std::cerr << "OGS_DIRICHLET_PERFORMANCE_GRID='" << grid
                      << "' is not a number, keeping " << n << '\n';
        }
    }
    // The triplet indices are cast to int, so n*n must stay within int.
    Eigen::Index const max_n =
        static_cast<Eigen::Index>(std::sqrt(std::numeric_limits<int>::max()));
    n = std::clamp<Eigen::Index>(n, 2, max_n);
    Eigen::Index const n_rows = n * n;
    int const repetitions = 5;

    std::mt19937 gen{42};
    std::uniform_real_distribution dist{1., 2.};

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(9 * n_rows);
    for (Eigen::Index i = 0; i < n; ++i)
    {
        for (Eigen::Index j = 0; j < n; ++j)
        {
            for (Eigen::Index di = -1; di <= 1; ++di)
            {
                for (Eigen::Index dj = -1; dj <= 1; ++dj)
                {
                    if (i + di < 0 || i + di >= n || j + dj < 0 || j + dj >= n)
                    {
                        continue;
                    }
                    triplets.emplace_back(
                        static_cast<int>(i * n + j),
                        static_cast<int>((i + di) * n + j + dj), dist(gen));
                }
            }
        }
    }
    MathLib::EigenMatrix A{n_rows};
    A.getRawMatrix().setFromTriplets(triplets.begin(), triplets.end());
    A.getRawMatrix().makeCompressed();

    MathLib::EigenVector rhs{n_rows};
    for (Eigen::Index row = 0; row < n_rows; ++row)
    {
        rhs[row] = dist(gen);
    }

    std::vector<Eigen::Index> known_x_ids;
    for (Eigen::Index i = 0; i < n; ++i)
    {
        for (Eigen::Index j = 0; j < n; ++j)
        {
            if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
            {
                known_x_ids.push_back(i * n + j);
            }
        }
    }
    std::vector<double> known_x_values(known_x_ids.size());
    for (auto& known_x : known_x_values)
    {
        known_x = dist(gen);
    }

    auto time = [&](auto&& apply)
    {
        double best = std::numeric_limits<double>::max();
        for (int r = 0; r < repetitions; ++r)
        {
            auto A_copy = A;
            auto rhs_copy = rhs;
            auto const start = std::chrono::steady_clock::now();
            apply(A_copy, rhs_copy);
            std::chrono::duration<double> const elapsed =
                std::chrono::steady_clock::now() - start;
            best = std::min(best, elapsed.count());
        }
        return best;
    };

    double const t_former = time(
        [&](MathLib::EigenMatrix& A_copy, MathLib::EigenVector& rhs_copy)
        {
            applyKnownSolutionCompleteByTransposition(
                A_copy.getRawMatrix(), rhs_copy.getRawVector(), known_x_ids,
                known_x_values);
        });
    MathLib::EigenVector x{n_rows};
    double const t_new = time(
        [&](MathLib::EigenMatrix& A_copy, MathLib::EigenVector& rhs_copy)
        {
            MathLib::applyKnownSolution(
                A_copy, rhs_copy, x, known_x_ids, known_x_values,
                MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE);
        });

    std::cout << "rows: " << n_rows
              << ", nonzeros: " << A.getRawMatrix().nonZeros()
              << ", Dirichlet DOFs: " << known_x_ids.size() << '\n'
              << "former (transposition): " << t_former << " s\n"
              << "single pass:            " << t_new << " s\n"
              << "speed-up:               " << t_former / t_new << '\n';

    // Same result as the former implementation.
    auto A_former = A;
    auto rhs_former = rhs;
    applyKnownSolutionCompleteByTransposition(A_former.getRawMatrix(),
                                              rhs_former.getRawVector(),
                                              known_x_ids, known_x_values);
    auto A_new = A;
    auto rhs_new = rhs;
    MathLib::applyKnownSolution(
        A_new, rhs_new, x, known_x_ids, known_x_values,
        MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE);
    EXPECT_TRUE(A_former.getRawMatrix().isApprox(A_new.getRawMatrix(), 0.0));
    EXPECT_TRUE(rhs_former.getRawVector() == rhs_new.getRawVector());
}
#endif
