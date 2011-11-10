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

namespace MathLib {

template<class T> class CRSMatrixOpenMP : public CRSMatrix<T, unsigned> {
public:
	CRSMatrixOpenMP(std::string const &fname, unsigned num_of_threads) :
			CRSMatrix<T, unsigned>(fname), _num_of_threads (num_of_threads)
	{}

	CRSMatrixOpenMP(unsigned n, unsigned *iA, unsigned *jA, T* A, unsigned num_of_threads) :
		CRSMatrix<T, unsigned>(n, iA, jA, A), _num_of_threads (num_of_threads)
	{}

	CRSMatrixOpenMP(unsigned n1) :
		CRSMatrix<T, unsigned>(n1), _num_of_threads (1)
	{}

	virtual ~CRSMatrixOpenMP()
	{}

	virtual void amux(T d, T const * const x, T *y) const
	{
		amuxCRSParallelOpenMP(d, MatrixBase::_n_rows, CRSMatrix<T,unsigned>::_row_ptr, CRSMatrix<T,unsigned>::_col_idx, CRSMatrix<T,unsigned>::_data, x, y, _num_of_threads);
	}

private:
	unsigned _num_of_threads;
};

} // end namespace MathLib

#endif /* CRSMATRIXOPENMP_H_ */
