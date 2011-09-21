#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <string>
#include <fstream>
#include <iostream>
#include "SparseMatrixBase.h"
#include "sparse.h"
#include "amuxCRS.h"

template<class T> class CRSMatrix : public SparseMatrixBase<T>
{
public:
	CRSMatrix(std::string const &fname) :
		SparseMatrixBase<T>(),
		_row_ptr(NULL), _col_idx(NULL), _data(NULL)
	{
		std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
		if (in) {
			CS_read(in, SparseMatrixBase<T>::_n_rows, _row_ptr, _col_idx, _data);
			SparseMatrixBase<T>::_n_cols = SparseMatrixBase<T>::_n_rows;
			in.close();
		} else {
			std::cout << "cannot open " << fname << std::endl;
		}
	}

	CRSMatrix(unsigned n, unsigned *iA, unsigned *jA, T* A) :
		SparseMatrixBase<T>(n,n), _row_ptr(iA), _col_idx(jA), _data(A)
	{}

	CRSMatrix(unsigned n1) :
		SparseMatrixBase<T>(n1, n1), _row_ptr(NULL), _col_idx(NULL), _data(NULL)
	{}

	virtual ~CRSMatrix()
	{
		delete [] _row_ptr;
		delete [] _col_idx;
		delete [] _data;
	}

	virtual void amux(T d, T const * const x, T *y) const
	{
		amuxCRS(d, SparseMatrixBase<T>::_n_rows, _row_ptr, _col_idx, _data, x, y);
	}

protected:
	unsigned *_row_ptr;
	unsigned *_col_idx;
	T* _data;
};

#endif

