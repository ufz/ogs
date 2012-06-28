/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file CRSMatrixPThreads.h
 *
 * Created on 2011-08-02 by Thomas Fischer
 */

#ifndef CRSMATRIXPTHREADS_H
#define CRSMATRIXPTHREADS_H

#include <string>

#include "SparseMatrixBase.h"
#include "sparse.h"
#include "CRSMatrix.h"
#include "amuxCRS.h"

namespace MathLib {

template<class T> class CRSMatrixPThreads : public CRSMatrix<T,unsigned>
{
public:
	CRSMatrixPThreads(std::string const &fname, unsigned num_of_threads) :
		CRSMatrix<T,unsigned>(fname), _num_of_threads (num_of_threads)
	{}

	CRSMatrixPThreads(unsigned n, unsigned *iA, unsigned *jA, T* A, unsigned num_of_threads) :
		CRSMatrix<T,unsigned>(n, iA, jA, A), _num_of_threads (num_of_threads)
	{}

	CRSMatrixPThreads(unsigned n1) :
		CRSMatrix<T,unsigned>(n1), _num_of_threads (1)
	{}

	virtual ~CRSMatrixPThreads()
	{}

	virtual void amux(T d, T const * const x, T *y) const
	{
		amuxCRSParallelPThreads(d, SparseMatrixBase<T, unsigned>::_n_rows,
						CRSMatrix<T, unsigned>::_row_ptr, CRSMatrix<T, unsigned>::_col_idx,
						CRSMatrix<T, unsigned>::_data, x, y, _num_of_threads);
	}

protected:
	unsigned _num_of_threads;
};

} // end namespace MathLib

#endif

