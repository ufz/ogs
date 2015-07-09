/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenTools.h"

#include <logog/include/logog.hpp>

#include "EigenMatrix.h"
#include "EigenVector.h"

namespace MathLib
{

void applyKnownSolution(EigenMatrix &A_, EigenVector &b_, const std::vector<std::size_t> &vec_knownX_id,
		const std::vector<double> &vec_knownX_x, double /*penalty_scaling*/)
{
    typedef EigenMatrix::RawMatrixType SpMat;
    auto &A = A_.getRawMatrix();
    auto &b = b_.getRawVector();
    const std::size_t n_rows = A.rows();
    for (std::size_t ix=0; ix<vec_knownX_id.size(); ix++)
    {
        int row_id = vec_knownX_id[ix];
        auto x = vec_knownX_x[ix];
        //A(k, j) = 0.
        for (SpMat::InnerIterator it(A,row_id); it; ++it)
            it.valueRef() = .0;
        //b_i -= A(i,k)*val, i!=k
        for (std::size_t i=0; i<n_rows; i++)
            for (SpMat::InnerIterator it(A,i); it; ++it)
            {
                if (it.col()!=row_id) continue;
                b[i] -= it.value()*x;
            }
        //b_k = val
        b[row_id] = x;
        //A(i, k) = 0., i!=k
        for (std::size_t i=0; i<n_rows; i++)
            for (SpMat::InnerIterator it(A,i); it; ++it)
            {
                if (it.col()!=row_id) continue;
                it.valueRef() = 0.0;
            }
        //A(k, k) = 1.0
        A.coeffRef(row_id, row_id) = 1.0; //=x
    }
}

} // MathLib



