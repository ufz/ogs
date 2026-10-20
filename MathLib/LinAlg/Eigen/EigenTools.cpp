// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "EigenTools.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <range/v3/view/zip.hpp>
#include <utility>

#include "EigenVector.h"

namespace
{
void setRowsToZeroOffDiagonal(
    MathLib::EigenMatrix::RawMatrixType& A_eigen,
    const std::vector<MathLib::EigenMatrix::IndexType>& row_ids)
{
    using SpMat = MathLib::EigenMatrix::RawMatrixType;

    // A_eigen(k, j) = 0.
    // set row to zero
    for (auto row_id : row_ids)
    {
        for (SpMat::InnerIterator it(A_eigen, row_id); it; ++it)
        {
            if (it.col() != decltype(it.col())(row_id))
            {
                it.valueRef() = 0.0;
            }
        }
    }
}

/// Symmetric elimination of the known solutions: the off-diagonal entries of
/// the known rows and columns are zeroed, the latter after their
/// contributions have been moved to the right-hand side, and the known rows
/// are scaled to the diagonal entry.
///
/// The columns are processed in a single pass over the row-major matrix with
/// a dense marker of the known DOFs. This avoids transposing the matrix twice,
/// i.e., two copies of the whole matrix. The sparsity pattern is kept, except
/// that missing diagonal entries of known rows are inserted.
///
/// The contributions of the known columns are subtracted from each row of
/// the right-hand side in the order of \c vec_knownX_id, so the result is
/// identical to the one of the former implementation via the transposed
/// matrix whenever its "b -= a*x" is fused to a single rounding. Where the
/// build leaves "b -= a*x" uncontracted, the former implementation rounds
/// twice per product while std::fma below rounds once, and the results
/// differ by at most one rounding per product. If an id occurs multiple
/// times, e.g., for a node with a Dirichlet BC that also lies in a
/// deactivated subdomain, the value of its first occurrence is moved to the
/// right-hand side, while the known row gets the value of its last
/// occurrence. The ids must be valid row indices of \c A_eigen and
/// \c vec_knownX_x must hold one value per id, otherwise OGS_FATAL below
/// throws std::runtime_error.
void applyKnownSolutionComplete(
    MathLib::EigenMatrix::RawMatrixType& A_eigen,
    MathLib::EigenVector::RawVectorType& b_eigen,
    const std::vector<MathLib::EigenMatrix::IndexType>& vec_knownX_id,
    const std::vector<double>& vec_knownX_x)
{
    using SpMat = MathLib::EigenMatrix::RawMatrixType;

    if (vec_knownX_id.empty())
    {
        return;
    }

    // Position of the first occurrence of each known id in vec_knownX_id,
    // sized for the columns as well, which the single pass also looks up.
    constexpr std::size_t not_known = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> first_position(
        std::max(A_eigen.rows(), A_eigen.cols()), not_known);
    for (std::size_t i = 0; i < vec_knownX_id.size(); ++i)
    {
        auto const id = vec_knownX_id[i];
        if (id < 0 || id >= A_eigen.rows())
        {
            OGS_FATAL("Known-solution id {} is not a row of the matrix.", id);
        }
        auto& position = first_position[id];
        if (position == not_known)
        {
            position = i;
        }
    }

    // Contributions of the known columns k to a row i: the value of
    // A_eigen(i,k) and the position of k in vec_knownX_id.
    std::vector<std::pair<std::size_t, double>> contributions;
    for (SpMat::Index row_id = 0; row_id < A_eigen.outerSize(); ++row_id)
    {
        bool const is_known_row = first_position[row_id] != not_known;
        contributions.clear();
        for (SpMat::InnerIterator it(A_eigen, row_id); it; ++it)
        {
            auto const col_id = it.col();
            if (col_id == row_id)
            {
                continue;
            }
            if (is_known_row)
            {
                // Known rows do not contribute to the right-hand side; they
                // are overwritten below.
                it.valueRef() = 0.0;
                continue;
            }
            auto const position = first_position[col_id];
            if (position != not_known)
            {
                contributions.emplace_back(position, it.value());
                it.valueRef() = 0.0;
            }
        }

        // b_i -= A_eigen(i,k)*x_k, i!=k, in the order of vec_knownX_id. The
        // product and the subtraction must round exactly as they do in the
        // test's reference applyKnownSolutionCompleteByTransposition. Both
        // implementations accumulate each product with std::fma, single
        // fused rounding, so both round identically whatever the build
        // does with "b -= a*x".
        std::sort(contributions.begin(), contributions.end(),
                  [](auto const& a, auto const& b)
                  { return a.first < b.first; });
        for (auto const& contribution : contributions)
        {
            b_eigen[row_id] = std::fma(-contribution.second,
                                       vec_knownX_x[contribution.first],
                                       b_eigen[row_id]);
        }
    }

    for (auto const& [row_id, x] :
         ranges::views::zip(vec_knownX_id, vec_knownX_x))
    {
        // Inserts the diagonal entry if missing, e.g., for rows of
        // deactivated subdomains, which are empty.
        auto& c = A_eigen.coeffRef(row_id, row_id);
        if (c != 0.0)
        {
            b_eigen[row_id] = x * c;
        }
        else
        {
            b_eigen[row_id] = x;
            c = 1.0;
        }
    }

    if (!A_eigen.isCompressed())
    {
        A_eigen.makeCompressed();
    }
}

void applyKnownSolutionIncomplete(
    MathLib::EigenMatrix::RawMatrixType& A_eigen,
    MathLib::EigenVector::RawVectorType& b_eigen,
    const std::vector<MathLib::EigenMatrix::IndexType>& vec_knownX_id,
    const std::vector<double>& vec_knownX_x)
{
    INFO("partial Dirichlet BC application (rhs only).");

    setRowsToZeroOffDiagonal(A_eigen, vec_knownX_id);

    Eigen::VectorXd x_known{A_eigen.rows()};
    x_known.setZero();

    for (auto const& [row_id, x] :
         ranges::views::zip(vec_knownX_id, vec_knownX_x))
    {
        x_known[row_id] = x;

        // set rhs to Dirichlet value
        auto const c = A_eigen.coeff(row_id, row_id);
        if (c != 0)
        {
            b_eigen[row_id] = x * c;

            // exclude diagonal (A_eigen(i,i)) from multiplication below
            // (A_eigen * x_known)
            A_eigen.coeffRef(row_id, row_id) = 0;
        }
        else
        {
            b_eigen[row_id] = x;
        }
    }

    b_eigen -= A_eigen * x_known;
}
}  // namespace

namespace MathLib
{
void applyKnownSolution(
    EigenMatrix& A, EigenVector& b, EigenVector& /*x*/,
    const std::vector<EigenMatrix::IndexType>& vec_knownX_id,
    const std::vector<double>& vec_knownX_x,
    DirichletBCApplicationMode const mode)
{
    if (vec_knownX_id.size() != vec_knownX_x.size())
    {
        OGS_FATAL(
            "Known-solution ids and values must pair up: {} ids, {} values.",
            vec_knownX_id.size(), vec_knownX_x.size());
    }

    using SpMat = EigenMatrix::RawMatrixType;
    static_assert(SpMat::IsRowMajor, "matrix is assumed to be row major!");

    auto& A_eigen = A.getRawMatrix();
    auto& b_eigen = b.getRawVector();

    using enum DirichletBCApplicationMode;
    switch (mode)
    {
        case COMPLETE_MATRIX_UPDATE:
            ::applyKnownSolutionComplete(A_eigen, b_eigen, vec_knownX_id,
                                         vec_knownX_x);
            return;
        case FAST_INCOMPLETE_MATRIX_UPDATE:
            ::applyKnownSolutionIncomplete(A_eigen, b_eigen, vec_knownX_id,
                                           vec_knownX_x);
            return;
    }

    OGS_FATAL("Unhandled DirichletBCApplicationMode with integer value {}.",
              std::to_underlying(mode));
}

}  // namespace MathLib
