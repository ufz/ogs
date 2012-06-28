/**
 * \file CRSSymMatrix.h
 *
 * Created on 2011-07-22 by Thomas Fischer
 */

#ifndef CRSSYMMATRIX_H_
#define CRSSYMMATRIX_H_

#include "CRSMatrix.h"

template<class T> class CRSSymMatrix : public CRSMatrix<T>
{
public:
	CRSSymMatrix(std::string const &fname)
	: CRSMatrix<T> (fname)
	{
		unsigned nnz (0);

		// count number of non-zeros in the upper triangular part
		for (unsigned i = 0; i < SparseMatrixBase<T>::_n_rows; i++) {
			unsigned idx = CRSMatrix<T>::_row_ptr[i+1];
			for (unsigned j = CRSMatrix<T>::_row_ptr[i]; j < idx; j++)
				if (CRSMatrix<T>::_col_idx[j] >= i)
					++nnz;
		}

		double *A_new (new double[nnz]);
		unsigned *jA_new (new unsigned[nnz]);
		unsigned *iA_new (new unsigned[SparseMatrixBase<T>::_n_rows+1]);

		iA_new[0] = nnz = 0;

		for (unsigned i = 0; i < SparseMatrixBase<T>::_n_rows; i++) {
			const unsigned idx (CRSMatrix<T>::_row_ptr[i+1]);
			for (unsigned j = CRSMatrix<T>::_row_ptr[i]; j < idx; j++) {
				if (CRSMatrix<T>::_col_idx[j] >= i) {
					A_new[nnz] = CRSMatrix<T>::_data[j];
					jA_new[nnz++] = CRSMatrix<T>::_col_idx[j];
				}
			}
			iA_new[i+1] = nnz;
		}

		std::swap(CRSMatrix<T>::_row_ptr, iA_new);
		std::swap(CRSMatrix<T>::_col_idx, jA_new);
		std::swap(CRSMatrix<T>::_data, A_new);

		delete[] iA_new;
		delete[] jA_new;
		delete[] A_new;
	}

	virtual ~CRSSymMatrix() {}

	void amux(T d, T const * const x, T *y) const
	{
		amuxCRSSym (d, SparseMatrixBase<T>::_n_rows, CRSMatrix<T>::_row_ptr, CRSMatrix<T>::_col_idx, CRSMatrix<T>::_data, x, y);
	}

};

#endif /* CRSSYMMATRIX_H_ */
