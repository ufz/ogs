/*
 * CRSMatrixOpenMP.h
 *
 *  Created on: Aug 8, 2011
 *      Author: TF
 */

#ifndef CRSMATRIXOPENMP_H_
#define CRSMATRIXOPENMP_H_

#include <string>

#include "CRSMatrix.h"
#include "amuxCRS.h"

template<class T> class CRSMatrixOpenMP : public CRSMatrix<T> {
public:
	CRSMatrixOpenMP(std::string const &fname, unsigned num_of_threads) :
			CRSMatrix<T>(fname), _num_of_threads (num_of_threads)
	{}

	CRSMatrixOpenMP(unsigned n, unsigned *iA, unsigned *jA, T* A, unsigned num_of_threads) :
		CRSMatrix<T>(n, iA, jA, A), _num_of_threads (num_of_threads)
	{}

	CRSMatrixOpenMP(unsigned n1) :
		CRSMatrix<T>(n1), _num_of_threads (1)
	{}

	virtual ~CRSMatrixOpenMP()
	{}

	virtual void amux(T d, T const * const x, T *y) const
	{
		amuxCRSParallelOpenMP(d, SparseMatrixBase<T>::_n_rows, CRSMatrix<T>::_row_ptr, CRSMatrix<T>::_col_idx, CRSMatrix<T>::_data, x, y, _num_of_threads);
	}

private:
	unsigned _num_of_threads;
};

#endif /* CRSMATRIXOPENMP_H_ */
