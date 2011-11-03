/*
 * CRSMatrix.h
 *
 *  Created on: Sep 20, 2011
 *      Author: TF
 */

#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include "SparseMatrixBase.h"
#include "sparse.h"
#include "amuxCRS.h"
#include "../Preconditioner/generateDiagPrecond.h"

namespace MathLib {

template<typename FP_TYPE, typename IDX_TYPE>
class CRSMatrix: public SparseMatrixBase<FP_TYPE, IDX_TYPE>
{
public:
	CRSMatrix(std::string const &fname) :
		SparseMatrixBase<FP_TYPE, IDX_TYPE>(),
		_row_ptr(NULL), _col_idx(NULL), _data(NULL)
	{
		std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
		if (in) {
			CS_read(in, SparseMatrixBase<FP_TYPE, IDX_TYPE>::_n_rows, _row_ptr, _col_idx, _data);
			SparseMatrixBase<FP_TYPE, IDX_TYPE>::_n_cols = SparseMatrixBase<FP_TYPE, IDX_TYPE>::_n_rows;
			in.close();
		} else {
			std::cout << "cannot open " << fname << std::endl;
		}
	}

	CRSMatrix(IDX_TYPE n, IDX_TYPE *iA, IDX_TYPE *jA, FP_TYPE* A) :
		SparseMatrixBase<FP_TYPE, IDX_TYPE>(n,n), _row_ptr(iA), _col_idx(jA), _data(A)
	{}

	CRSMatrix(IDX_TYPE n1) :
		SparseMatrixBase<FP_TYPE, IDX_TYPE>(n1, n1), _row_ptr(NULL), _col_idx(NULL), _data(NULL)
	{}

	virtual ~CRSMatrix()
	{
		delete [] _row_ptr;
		delete [] _col_idx;
		delete [] _data;
	}

	virtual void amux(FP_TYPE d, FP_TYPE const * const x, FP_TYPE *y) const
	{
		amuxCRS(d, MatrixBase::_n_rows, _row_ptr, _col_idx, _data, x, y);
	}

    virtual void precondApply(FP_TYPE* x) const
    {}

    /**
     * get the number of non-zero entries
     * @return number of non-zero entries
     */
    IDX_TYPE getNNZ() const { return _row_ptr[MatrixBase::_n_rows]; }

    /**
     * This is the constant access operator to a non-zero matrix entry.
     * Precondition: the entries have to be in the sparsity pattern!
     * @param row the row number
     * @param col the column number
     * @return The corresponding matrix entry.
     */
    FP_TYPE const& operator() (IDX_TYPE row, IDX_TYPE col) const
    {
    	assert(0 <= row && row < MatrixBase::_n_rows);

    	// linear search - for matrices with many entries per row binary search is much faster
    	const IDX_TYPE idx_end (_row_ptr[row+1]);
    	IDX_TYPE j(_row_ptr[row]), k;

    	while (j<idx_end && (k=_col_idx[j]) <= col) {
    		if (k == col) {
    			return _data[j];
    		}
    		j++;
    	}
    }

    /**
	 * This is the non constant access operator to a non-zero matrix entry.
	 * Precondition: the entries have to be in the sparsity pattern!
	 * @param row the row number
	 * @param col the column number
	 * @return The corresponding matrix entry that can be modified by the user.
	 */
	FP_TYPE & operator() (IDX_TYPE row, IDX_TYPE col)
	{
		assert(0 <= row && row < MatrixBase::_n_rows);

		// linear search - for matrices with many entries per row binary search is much faster
		const IDX_TYPE idx_end (_row_ptr[row+1]);
		IDX_TYPE j(_row_ptr[row]), k;

		while (j<idx_end && (k=_col_idx[j]) <= col) {
			if (k == col) {
				return _data[j];
			}
			j++;
		}
	}

protected:
	IDX_TYPE *_row_ptr;
	IDX_TYPE *_col_idx;
	FP_TYPE* _data;
};

} // end namespace MathLib

#endif

