/**
 * \file
 * \author Thomas Fischer
 * \date   2011-08-02
 * \brief  Definition of the CRSMatrixPThreads class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
		CRSMatrix<T,unsigned>(fname), _n_threads (num_of_threads),
		_workload_intervals(new unsigned[num_of_threads+1])
	{
		calcWorkload();
	}

	CRSMatrixPThreads(unsigned n, unsigned *iA, unsigned *jA, T* A, unsigned num_of_threads) :
		CRSMatrix<T,unsigned>(n, iA, jA, A), _n_threads (num_of_threads),
		_workload_intervals(new unsigned[num_of_threads+1])
	{
		calcWorkload();
	}

	CRSMatrixPThreads(unsigned n1) :
		CRSMatrix<T,unsigned>(n1), _n_threads (1),
		_workload_intervals(new unsigned[_n_threads+1])
	{
		calcWorkload();
	}

	virtual ~CRSMatrixPThreads()
	{
		delete [] _workload_intervals;
	}

	virtual void amux(T d, T const * const x, T *y) const
	{
		amuxCRSParallelPThreads(d, SparseMatrixBase<T, unsigned>::_n_rows,
						CRSMatrix<T, unsigned>::_row_ptr, CRSMatrix<T, unsigned>::_col_idx,
						CRSMatrix<T, unsigned>::_data, x, y, _n_threads, _workload_intervals);
	}

protected:
	void calcWorkload()
	{
		_workload_intervals[0] = 0;
		_workload_intervals[_n_threads] = SparseMatrixBase<T, unsigned>::_n_rows;

		const unsigned work_per_core (this->getNNZ()/_n_threads);
		for (unsigned k(1); k<_n_threads; k++) {
			unsigned upper_bound_kth_core(k * work_per_core);
			// search in _row_ptr array for the appropriate index
			unsigned beg (_workload_intervals[k-1]);
			unsigned end (_workload_intervals[_n_threads]);
			bool found (false);
			while (beg < end && !found) {
				unsigned m ((end+beg)/2);

				if (upper_bound_kth_core == this->_row_ptr[m]) {
					_workload_intervals[k] = m;
					found = true;
				} else {
					if (upper_bound_kth_core < this->_row_ptr[m]) {
						end = m;
					} else {
						beg = m+1;
					}
				}
			}
			if (!found)
				_workload_intervals[k] = beg;
		}

		for (unsigned k(0); k<_n_threads; k++) {
			std::cout << "proc " << k << ": [" << _workload_intervals[k] << "," << _workload_intervals[k+1] << ") - "
				<< _workload_intervals[k+1] - _workload_intervals[k] << " rows and "
				<< this->_row_ptr[_workload_intervals[k+1]] - this->_row_ptr[_workload_intervals[k]] << " entries" << std::endl;
		}
	}

	const unsigned _n_threads;
	unsigned *_workload_intervals;
};

} // end namespace MathLib

#endif

