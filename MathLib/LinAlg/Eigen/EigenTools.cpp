/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenTools.h"

#include "EigenVector.h"

namespace MathLib
{
void applyKnownSolution(
    EigenMatrix& A, EigenVector& b, EigenVector& /*x*/,
    const std::vector<EigenMatrix::IndexType>& vec_knownX_id,
    const std::vector<double>& vec_knownX_x, double /*penalty_scaling*/)
{
    using SpMat = EigenMatrix::RawMatrixType;
    static_assert(SpMat::IsRowMajor, "matrix is assumed to be row major!");

    auto& A_eigen = A.getRawMatrix();
    auto& b_eigen = b.getRawVector();

    // A_eigen(k, j) = 0.
    // set row to zero
    for (auto row_id : vec_knownX_id)
    {
        for (SpMat::InnerIterator it(A_eigen, row_id); it; ++it)
        {
            if (it.col() != decltype(it.col())(row_id))
            {
                it.valueRef() = 0.0;
            }
        }
    }

    SpMat AT = A_eigen.transpose();

    // Reserve space for at least one value (on the diagonal). For deactivated
    // subdomains some rows and columns might end empty (in the A matrix) and so
    // no space is reserved for the Dirichlet conditions in the transposed
    // matrix. Then the coeffRef call will do costly reallocations.
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

}  // namespace MathLib
