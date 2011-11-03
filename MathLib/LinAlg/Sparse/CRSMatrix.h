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

    IDX_TYPE getNNZ() const { return _row_ptr[MatrixBase::_n_rows]; }

protected:
	IDX_TYPE *_row_ptr;
	IDX_TYPE *_col_idx;
	FP_TYPE* _data;
};

} // end namespace MathLib

#endif

