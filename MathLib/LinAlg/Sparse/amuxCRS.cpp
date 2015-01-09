/**
 * \file
 * \author Thomas Fischer
 * \date   2011-09-20
 * \brief  Implementation of amuxCRS functions.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "amuxCRS.h"
#include <cstddef>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_PTHREADS
#include <pthread.h>
#endif

namespace MathLib {

struct MatMultThreadParam {
	MatMultThreadParam (double scalar_factor, unsigned beg_row, unsigned end_row,
		unsigned const * const iA, unsigned const * const jA,
	        double const * const A, double const * const x, double* y) :
		_a (scalar_factor), _beg_row(beg_row), _end_row(end_row),
		_row_ptr(iA), _col_idx(jA), _data(A), _x(x), _y(y)
	{}

	double _a;
	unsigned _beg_row;
	unsigned _end_row;
	unsigned const * const _row_ptr;
	unsigned const * const _col_idx;
	double const * const _data;
	double const * const _x;
	double * _y;
};

extern "C" {
void* amuxCRSpthread (void* ptr)
{
	MatMultThreadParam *thread_param(static_cast<MatMultThreadParam*>(ptr));
	const double a(thread_param->_a);
	const unsigned beg_row(thread_param->_beg_row);
	const unsigned end_row(thread_param->_end_row);
	unsigned const * const iA(thread_param->_row_ptr);
	unsigned const * const jA(thread_param->_col_idx);
	double const * const A(thread_param->_data);
	double const * const x(thread_param->_x);
	double* y(thread_param->_y);

	for (unsigned i(beg_row); i<end_row; i++) {
		y[i] = A[iA[i]] * x[jA[iA[i]]];
		const unsigned end (iA[i+1]);
		for (unsigned j(iA[i]+1); j<end; j++) {
			y[i] += A[j] * x[jA[j]];
		}
		y[i] *= a;
	}
	return NULL;
}
} // end extern "C"

void amuxCRSParallelPThreads (double a,
	unsigned n, unsigned const * const iA, unsigned const * const jA,
	double const * const A, double const * const x, double* y,
	unsigned num_of_pthreads)
{
#ifdef HAVE_PTHREADS
	// fill thread data objects
	MatMultThreadParam** thread_param_array (new MatMultThreadParam*[num_of_pthreads]);
	double step_size (static_cast<double>(n)/(num_of_pthreads));
	for (unsigned k(0); k<num_of_pthreads; k++) {
		const unsigned beg (static_cast<unsigned>(k*step_size));
		const unsigned end (static_cast<unsigned>((k+1)*step_size));
		thread_param_array[k] = new MatMultThreadParam (a, beg, end, iA, jA, A, x, y);
	}

	// allocate thread_array and return value array
	pthread_t *thread_array (new pthread_t[num_of_pthreads]);
	int *ret_vals (new int[num_of_pthreads]);

	// create threads
	for (unsigned k(0); k<num_of_pthreads; k++) {
		ret_vals[k] = pthread_create( &(thread_array[k]), NULL, amuxCRSpthread, thread_param_array[k]);
	}

	// join threads
	for (unsigned k(0); k<num_of_pthreads; k++) {
		pthread_join (thread_array[k], NULL);
	}

	delete [] ret_vals;
	for (unsigned k(0); k<num_of_pthreads; k++)
		delete thread_param_array[k];
	delete [] thread_param_array;
	delete [] thread_array;
#else
	(void)num_of_pthreads;
	amuxCRS (a, n, iA, jA, A, x, y);
#endif
}

void amuxCRSParallelPThreads (double a,
	unsigned n, unsigned const * const iA, unsigned const * const jA,
	double const * const A, double const * const x, double* y,
	unsigned num_of_pthreads, unsigned const*const workload_intervals)
{
    (void) n;   // Unused if HAVE_PTHREADS is not defined.

#ifdef HAVE_PTHREADS
	// fill thread data objects
	MatMultThreadParam** thread_param_array (new MatMultThreadParam*[num_of_pthreads]);
	for (unsigned k(0); k<num_of_pthreads; k++) {
		thread_param_array[k] = new MatMultThreadParam (a, workload_intervals[k], workload_intervals[k+1], iA, jA, A, x, y);
	}

	// allocate thread_array and return value array
	pthread_t *thread_array (new pthread_t[num_of_pthreads]);
	int *ret_vals (new int[num_of_pthreads]);

	// create threads
	for (unsigned k(0); k<num_of_pthreads; k++) {
		ret_vals[k] = pthread_create( &(thread_array[k]), NULL, amuxCRSpthread, thread_param_array[k]);
	}

	// join threads
	for (unsigned k(0); k<num_of_pthreads; k++) {
		pthread_join (thread_array[k], NULL);
	}

	delete [] ret_vals;
	for (unsigned k(0); k<num_of_pthreads; k++)
		delete thread_param_array[k];
	delete [] thread_param_array;
	delete [] thread_array;
#else
	(void)num_of_pthreads;
	amuxCRS (a, n, iA, jA, A, x, y);
#endif
}


void amuxCRSSym (double a,
	unsigned n, unsigned const * const iA, unsigned const * const jA,
        double const * const A, double const * const x, double* y)
{
	for (unsigned i(0); i<n; i++) {
			y[i] = 0.0;
	}

	for (unsigned i(0); i<n; i++) {
		unsigned j (iA[i]);
		// handle diagonal
		if (jA[j] == i) {
			y[i] += A[j] * x[jA[j]];
			j++;
		}
		const unsigned end (iA[i+1]);
		for (; j<end; j++) {
				y[i] += A[j] * x[jA[j]];
				y[jA[j]] += A[j] * x[i];
		}
	}

	for (unsigned i(0); i<n; i++) {
		y[i] *= a;
	}
}

} // end namespace MathLib
