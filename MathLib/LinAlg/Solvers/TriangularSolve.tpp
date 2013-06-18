/**
 * \file
 * \author Thomas Fischer
 * \date   2011-05-05
 * \brief  Implementation of triangular solver functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace MathLib {

template <typename FP_T>
void forwardSolve (const DenseMatrix <FP_T> &L, FP_T* b)
{
	size_t m (L.getNRows());
	FP_T t;

	for (size_t r=0; r<m; r++) {
		t = 0.0;
		for (size_t c=0; c<r; c++) {
			t += L(r,c)*b[c];
		}
		b[r] = b[r]-t;
	}
}

template <typename FP_T>
void backwardSolve (const DenseMatrix <FP_T> &mat, FP_T* b)
{
	FP_T t;
	size_t m (mat.getNRows()), n(mat.getNCols());
	for (int r=m-1; r>=0; r--) {
		t = 0.0;
		for (size_t c=r+1; c<n; c++) {
			t += mat(r,c)*b[c];
		}
		b[r] = (b[r]-t) / mat(r,r);
	}
}

template <typename FP_T>
void backwardSolve ( DenseMatrix<FP_T> const& mat, FP_T* x, FP_T* b)
{
	size_t n_cols (mat.getNCols());
	for (int r = (n_cols - 1); r >= 0; r--) {
		FP_T t = 0.0;

		for (size_t c = r+1; c < n_cols; c++) {
			t += mat(r,c) * b[c];
		}
		x[r] = (b[r] - t) / mat(r, r);
	}
}


} // end namespace MathLib
