/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenTools.h"

#include <logog/include/logog.hpp>

#include "EigenVector.h"

namespace MathLib
{

void applyKnownSolution(EigenMatrix &A_, EigenVector &b_, EigenVector &/*x*/,
        const std::vector<EigenMatrix::IndexType> &vec_knownX_id,
        const std::vector<double> &vec_knownX_x, double /*penalty_scaling*/)
{
    using SpMat = EigenMatrix::RawMatrixType;
    static_assert(SpMat::IsRowMajor, "matrix is assumed to be row major!");

    auto &A = A_.getRawMatrix();
    auto &b = b_.getRawVector();

    // A(k, j) = 0.
    // set row to zero
    for (auto row_id : vec_knownX_id)
        for (SpMat::InnerIterator it(A,row_id); it; ++it) {
            if (it.col() != decltype(it.col())(row_id)) it.valueRef() = 0.0;
        }

    SpMat AT = A.transpose();

    for (std::size_t ix=0; ix<vec_knownX_id.size(); ix++)
    {
        SpMat::Index const row_id = vec_knownX_id[ix];
        auto const x = vec_knownX_x[ix];

        // b_i -= A(i,k)*val, i!=k
        // set column to zero, subtract from rhs
        for (SpMat::InnerIterator it(AT, row_id); it; ++it)
        {
            if (it.col() == row_id) continue;

            b[it.col()] -= it.value()*x;
            it.valueRef() = 0.0;
        }

        auto& c = AT.coeffRef(row_id, row_id);
        if (c != 0.0) {
            b[row_id] = x * c;
        } else {
            b[row_id] = x;
            c = 1.0;
        }
    }

    A = AT.transpose();
}

} // MathLib



