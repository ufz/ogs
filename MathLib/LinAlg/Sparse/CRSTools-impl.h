/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>
#include "logog/include/logog.hpp"

#include "CRSTools.h"

#include "CRSMatrix.h"

namespace MathLib
{

template <typename VEC_T, typename FP_TYPE>
void applyKnownSolution(CRSMatrix<FP_TYPE, typename VEC_T::IndexType>*& mat,
	VEC_T &rhs, std::vector<typename VEC_T::IndexType> const& rows,
	std::vector<FP_TYPE> const& vals)
{
	std::unique_ptr<MathLib::CRSMatrix<FP_TYPE, typename VEC_T::IndexType>> mat_t(mat->getTranspose());

	// row pointer of the transposed matrix
	typename VEC_T::IndexType const*const iAt(mat_t->getRowPtrArray());
	// column indices of the transposed matrix
	typename VEC_T::IndexType const*const jAt(mat_t->getColIdxArray());
	// entries of the transposed matrix
	double * At(const_cast<double *>(mat_t->getEntryArray()));

	// b_i -= A(i,k)*val, i!=k => b_i -= A(k,i)^T * val
	// set A^T(k,i) = 0, i!=k (i.e. set column entries of original matrix to
	// zero)
	for (std::size_t r(0); r<rows.size(); ++r) {
		auto const row = rows[r];
		auto const val = vals[r];
		for (typename VEC_T::IndexType j(iAt[row]); j<iAt[row+1]; ++j) {
			if (jAt[j] == row) // skip diagonal entry
				continue;
			rhs.add(jAt[j], -At[j] * val);
			At[j] = 0.0;
		}
	}

	delete mat;
	mat = mat_t->getTranspose();
	typename VEC_T::IndexType const*const iA(mat->getRowPtrArray()); // row ptrs
	typename VEC_T::IndexType const*const jA(mat->getColIdxArray()); // col idx
	double * entries(const_cast<double*>(mat->getEntryArray()));

	// set row entries, except the diagonal entry, to zero
	for (std::size_t r(0); r<rows.size(); ++r) {
		auto const row = rows[r];
		for (typename VEC_T::IndexType j = iA[row]; j < iA[row + 1]; ++j)
		{
			if (jA[j] == row) {
				entries[j] = 1.0; // A(row,row) = 1.0
				// rhs[row] = A(row,row) * vals[r]
				rhs.set(row, vals[r]);
			} else
				entries[j] = 0.0;
		}
	}
}

} // MathLib

