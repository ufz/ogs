// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "EigenTools.h"

#include <range/v3/view/zip.hpp>

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

void applyKnownSolutionComplete(
    MathLib::EigenMatrix::RawMatrixType& A_eigen,
    MathLib::EigenVector::RawVectorType& b_eigen,
    const std::vector<MathLib::EigenMatrix::IndexType>& vec_knownX_id,
    const std::vector<double>& vec_knownX_x)
{
    setRowsToZeroOffDiagonal(A_eigen, vec_knownX_id);

    using SpMat = MathLib::EigenMatrix::RawMatrixType;
    SpMat AT = A_eigen.transpose();

    // Reserve space for at least one value (on the diagonal). For
    // deactivated subdomains some rows and columns might end empty (in the
    // A matrix) and so no space is reserved for the Dirichlet conditions in
    // the transposed matrix. Then the coeffRef call will do costly
    // reallocations.
    AT.reserve(Eigen::VectorXi::Constant(A_eigen.rows(), 1));

    for (std::size_t ix = 0; ix < vec_knownX_id.size(); ix++)
    {
        SpMat::Index const row_id = vec_knownX_id[ix];
        auto const x = vec_knownX_x[ix];

        // b_i -= A_eigen(i,k)*val, i!=k
        // set column to zero, subtract from rhs
        for (SpMat::InnerIterator it(AT, row_id); it; ++it)
        {
            if (it.col() == row_id)
            {
                continue;
            }

            b_eigen[it.col()] -= it.value() * x;
            it.valueRef() = 0.0;
        }

        auto& c = AT.coeffRef(row_id, row_id);
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

    A_eigen = AT.transpose();
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
