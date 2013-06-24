/**
 * \file
 * \author Thomas Fischer
 * \date   2011-05-06
 * \brief  Implementation of the GaussAlgorithm class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <algorithm>
#include "GaussAlgorithm.h"

namespace MathLib {

template <typename MAT_T, typename VEC_T>
GaussAlgorithm<MAT_T, VEC_T>::GaussAlgorithm(MAT_T &A,
		boost::property_tree::ptree const* const) :
		_mat(A), _n(_mat.getNRows()), _perm(new IDX_T[_n])
{
	IDX_T k, i, j, nr (_mat.getNRows()), nc(_mat.getNCols());
	FP_T l;

	for (k=0; k<nc; k++) {
		// search pivot
		FP_T t = std::abs(_mat(k, k));
		_perm[k] = k;
		for (i=k+1; i<nr; i++) {
			if (std::abs(_mat(i,k)) > t) {
				t = _mat(i,k);
				_perm[k] = i;
			}
		}

		// exchange rows
		if (_perm[k] != k) {
			for (j=0; j<nc; j++)
				std::swap (_mat(_perm[k],j), _mat(k,j));
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

template <typename MAT_T, typename VEC_T>
GaussAlgorithm<MAT_T, VEC_T>::~GaussAlgorithm()
{
	delete [] _perm;
}

template <typename MAT_T, typename VEC_T>
template <typename V>
void GaussAlgorithm<MAT_T, VEC_T>::solve (V & b) const
{
	permuteRHS (b);
	forwardSolve (_mat, b); // L z = b, b will be overwritten by z
	backwardSolve (_mat, b); // U x = z, b (z) will be overwritten by x
}

template <typename MAT_T, typename VEC_T>
void GaussAlgorithm<MAT_T, VEC_T>:: solve(FP_T const* & b) const
{
	permuteRHS (b);
	forwardSolve (_mat, b); // L z = b, b will be overwritten by z
	backwardSolve (_mat, b); // U x = z, b (z) will be overwritten by x
}

template <typename MAT_T, typename VEC_T>
void GaussAlgorithm<MAT_T, VEC_T>::solve (FP_T* & b) const
{
	permuteRHS (b);
	forwardSolve (_mat, b); // L z = b, b will be overwritten by z
	backwardSolve (_mat, b); // U x = z, b (z) will be overwritten by x
}

template <typename MAT_T, typename VEC_T>
void GaussAlgorithm<MAT_T, VEC_T>::solve (VEC_T & x, VEC_T const& b) const
{
	for (std::size_t k(0); k<_mat.getNRows(); k++)
		x[k] = b[k];
	solve(x);
}

template <typename MAT_T, typename VEC_T>
template <typename V>
void GaussAlgorithm<MAT_T, VEC_T>::permuteRHS (V & b) const
{
	for (IDX_T i=0; i<_n; i++) {
		if (_perm[i] != i)
			std::swap(b[i], b[_perm[i]]);
	}
}

template <typename MAT_T, typename VEC_T>
void GaussAlgorithm<MAT_T, VEC_T>::permuteRHS (VEC_T& b) const
{
	for (IDX_T i=0; i<_n; i++) {
		if (_perm[i] != i)
			std::swap(b[i], b[_perm[i]]);
	}
}

} // end namespace MathLib
