/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file GaussAlgorithm.cpp
 *
 * Created on 2011-05-06 by Thomas Fischer
 */

#include <cmath>
#include "GaussAlgorithm.h"
#include "swap.h"

namespace MathLib {

GaussAlgorithm::GaussAlgorithm (Matrix <double> &A) :
	_mat (A), _n(_mat.getNRows()), _perm (new size_t [_n])
{
	size_t k, i, j, nr (_mat.getNRows()), nc(_mat.getNCols());
	double l;

	for (k=0; k<nc; k++) {
		// search pivot
		double t = fabs(_mat(k, k));
		_perm[k] = k;
		for (i=k+1; i<nr; i++) {
			if (fabs(_mat(i,k)) > t) {
				t = _mat(i,k);
				_perm[k] = i;
			}
		}

		// exchange rows
		if (_perm[k] != k) {
			for (j=0; j<nc; j++) BaseLib::swap (_mat(_perm[k],j), _mat(k,j));
		}

		// eliminate
		for (i=k+1; i<nr; i++) {
			l=_mat(i,k)/_mat(k,k);
			for (j=k; j<nc; j++) {
				_mat(i,j) -= _mat(k,j) * l;
			}
			_mat(i,k) = l;
		}
	}
}

GaussAlgorithm::~GaussAlgorithm()
{
	delete [] _perm;
}

void GaussAlgorithm::execute (double *b) const
{
	permuteRHS (b);
	forwardSolve (_mat, b); // L z = b, b will be overwritten by z
	backwardSolve (_mat, b); // U x = z, b (z) will be overwritten by x
}

void GaussAlgorithm::permuteRHS (double* b) const
{
	for (size_t i=0; i<_n; i++) {
		if (_perm[i] != i) BaseLib::swap(b[i], b[_perm[i]]);
	}
}

} // end namespace MathLib
