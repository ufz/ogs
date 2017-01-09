/**
 * \file
 * \author Thomas Fischer
 * \date   2011-05-05
 * \brief  Implementation of triangular solver functions.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace MathLib {

template <typename FP_T, typename VEC_T>
void forwardSolve (const DenseMatrix <FP_T> &L, VEC_T& b)
{
    typedef typename DenseMatrix<FP_T>::IDX_T IDX_T;
    IDX_T m (L.getNumberOfRows());
    FP_T t;

    for (IDX_T r=0; r<m; r++) {
        t = 0.0;
        for (IDX_T c=0; c<r; c++) {
            t += L(r,c)*b[c];
        }
        b[r] = b[r]-t;
    }
}

template <typename FP_T, typename VEC_T>
void backwardSolve (const DenseMatrix <FP_T> &mat, VEC_T& b)
{
    FP_T t;
    typedef typename DenseMatrix<FP_T>::IDX_T IDX_T;
    IDX_T m (mat.getNumberOfRows()), n(mat.getNumberOfColumns());
    for (int r=m-1; r>=0; r--) {
        t = 0.0;
        for (IDX_T c=r+1; c<n; c++) {
            t += mat(r,c)*b[c];
        }
        b[r] = (b[r]-t) / mat(r,r);
    }
}

template <typename FP_T, typename VEC_T>
void backwardSolve ( DenseMatrix<FP_T> const& mat, VEC_T& x, VEC_T const& b)
{
    typedef typename DenseMatrix<FP_T>::IDX_T IDX_T;
    IDX_T n_cols (mat.getNumberOfColumns());
    for (int r = (n_cols - 1); r >= 0; r--) {
        FP_T t = 0.0;

        for (IDX_T c = r+1; c < n_cols; c++) {
            t += mat(r,c) * b[c];
        }
        x[r] = (b[r] - t) / mat(r, r);
    }
}


} // end namespace MathLib
