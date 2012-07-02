/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TriangularSolve.cpp
 *
 * Created on 2011-05-05 by Thomas Fischer
 */

#include "../Dense/Matrix.h"

namespace MathLib {

void forwardSolve (const Matrix <double> &L, double* b)
{
	size_t m (L.getNRows());
	double t;

	for (size_t r=0; r<m; r++) {
		t = 0.0;
		for (size_t c=0; c<r; c++) {
			t += L(r,c)*b[c];
		}
		b[r] = b[r]-t;
	}
}

void backwardSolve (const Matrix <double> &mat, double* b)
{
	double t;
	size_t m (mat.getNRows()), n(mat.getNCols());
	for (int r=m-1; r>=0; r--) {
		t = 0.0;
		for (size_t c=r+1; c<n; c++) {
			t += mat(r,c)*b[c];
		}
		b[r] = (b[r]-t) / mat(r,r);
	}
}

void backwardSolve ( Matrix<double> const& mat, double* x, double* b)
{
	size_t n_cols (mat.getNCols());
	for (int r = (n_cols - 1); r >= 0; r--) {
		double t = 0.0;

		for (size_t c = r+1; c < n_cols; c++) {
			t += mat(r,c) * b[c];
		}
		x[r] = (b[r] - t) / mat(r, r);
	}
}


} // end namespace MathLib
