/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CRSMatrixReordered.cpp
 *
 * Created on 2012-01-03 by Thomas Fischer
 */

// BaseLib
#include "quicksort.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/CRSMatrixReordered.h"

namespace MathLib {

CRSMatrixReordered::CRSMatrixReordered(std::string const &fname) :
	CRSMatrix<double,unsigned>(fname)
{}

CRSMatrixReordered::CRSMatrixReordered(unsigned n, unsigned *iA, unsigned *jA, double* A) :
	CRSMatrix<double, unsigned> (n, iA, jA, A)
{}

CRSMatrixReordered::~CRSMatrixReordered()
{}

void CRSMatrixReordered::reorderMatrix(unsigned const*const op_perm, unsigned const*const po_perm)
{
	unsigned i; // row and col idx in permuted matrix
	unsigned j, idx; // pointer in jA

	const unsigned size(getNRows());

	unsigned *pos(new unsigned[size + 1]);
	for (i = 0; i < size; i++) {
		const unsigned original_row(op_perm[i]);
		pos[i] = _row_ptr[original_row+1] - _row_ptr[original_row];
	}
	pos[size] = 0;

	unsigned *iAn(new unsigned[size + 1]);
	iAn[0] = 0;
	for (i = 0; i < size; i++)
		iAn[i + 1] = iAn[i] + pos[i];
	for (i = 0; i < size; i++)
		pos[i] = iAn[i];

	unsigned *jAn(new unsigned[iAn[size]]);
	double *An(new double[iAn[size]]);
	for (i = 0; i < size; i++) {
		const unsigned original_row(op_perm[i]);
		idx = _row_ptr[original_row+1];
		for (j = _row_ptr[original_row]; j < idx; j++) {
			jAn[pos[i]] = po_perm[_col_idx[j]];
			An[pos[i]++] = _data[j];
		}
	}

	delete[] pos;
	for (i = 0; i < size; ++i)
		BaseLib::quicksort(jAn, static_cast<size_t>(iAn[i]), static_cast<size_t>(iAn[i + 1]), An);

	BaseLib::swap(iAn, _row_ptr);
	BaseLib::swap(jAn, _col_idx);
	BaseLib::swap(An, _data);

	delete [] iAn;
	delete [] jAn;
	delete [] An;
}

} // end namespace MathLib
