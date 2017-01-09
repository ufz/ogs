/**
 * \file
 * \author Thomas Fischer
 * \date   2011-05-06
 * \brief  Implementation of the GaussAlgorithm class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <algorithm>
#include "GaussAlgorithm.h"

namespace MathLib
{

template <typename MAT_T, typename VEC_T>
void GaussAlgorithm<MAT_T, VEC_T>::performLU(MAT_T& A)
{
    IDX_T const nr(A.getNumberOfRows());
    IDX_T const nc(A.getNumberOfColumns());

    for (IDX_T k=0; k<nc; k++) {
        // search pivot
        FP_T t = std::abs(A(k, k));
        _perm[k] = k;
        for (IDX_T i=k+1; i<nr; i++) {
            FP_T const s = std::abs(A(i,k));
            if (s > t) {
                t = s;
                _perm[k] = i;
            }
        }

        // exchange rows
        if (_perm[k] != k) {
            for (IDX_T j=0; j<nc; j++)
                std::swap (A(_perm[k],j), A(k,j));
        }

        // eliminate
        for (IDX_T i=k+1; i<nr; i++) {
            FP_T const l = A(i,k)/A(k,k);
            for (IDX_T j=k; j<nc; j++) {
                A(i,j) -= A(k,j) * l;
            }
            A(i,k) = l;
        }
    }
}

template <typename MAT_T, typename VEC_T>
template <typename V>
void GaussAlgorithm<MAT_T, VEC_T>::
solve (MAT_T& A, V& b, bool decompose)
{
    _perm.resize(A.getNumberOfRows());

    if (decompose)
        performLU(A);
    permuteRHS (b);
    forwardSolve (A, b); // L z = b, b will be overwritten by z
    backwardSolve (A, b); // U x = z, b (z) will be overwritten by x
}

template <typename MAT_T, typename VEC_T>
void GaussAlgorithm<MAT_T, VEC_T>::
solve (MAT_T& A, FP_T* & b, bool decompose)
{
    _perm.resize(A.getNumberOfRows());

    if (decompose)
        performLU(A);
    permuteRHS (b);
    forwardSolve (A, b); // L z = b, b will be overwritten by z
    backwardSolve (A, b); // U x = z, b (z) will be overwritten by x
}

template <typename MAT_T, typename VEC_T>
void GaussAlgorithm<MAT_T, VEC_T>::solve (
        MAT_T& A, VEC_T const& b, VEC_T & x,
        bool decompose)
{
    for (std::size_t k(0); k<A.getNumberOfRows(); k++)
        x[k] = b[k];
    solve(A, x, decompose);
}

template <typename MAT_T, typename VEC_T>
template <typename V>
void GaussAlgorithm<MAT_T, VEC_T>::permuteRHS (V & b) const
{
    for (IDX_T i=0; i<_perm.size(); i++) {
        if (_perm[i] != i)
            std::swap(b[i], b[_perm[i]]);
    }
}

template <typename MAT_T, typename VEC_T>
void GaussAlgorithm<MAT_T, VEC_T>::permuteRHS (VEC_T& b) const
{
    for (IDX_T i=0; i<_perm.size(); i++) {
        if (_perm[i] != i)
            std::swap(b[i], b[_perm[i]]);
    }
}

} // end namespace MathLib
