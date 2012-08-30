/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file amuxCRS.h
 *
 * Created on 2011-09-20 by Thomas Fischer
 */

#ifndef AMUXCRS_H
#define AMUXCRS_H

namespace MathLib {

template<typename FP_TYPE, typename IDX_TYPE>
void amuxCRS(FP_TYPE a, IDX_TYPE n, IDX_TYPE const * const iA, IDX_TYPE const * const jA,
				FP_TYPE const * const A, FP_TYPE const * const x, FP_TYPE* y)
{
	for (IDX_TYPE i(0); i < n; i++) {
		const IDX_TYPE end(iA[i + 1]);
		y[i] = A[iA[i]] * x[jA[iA[i]]];
		for (IDX_TYPE j(iA[i]+1); j < end; j++) {
			y[i] += A[j] * x[jA[j]];
		}
		y[i] *= a;
	}
}

void amuxCRSParallelPThreads (double a,
	unsigned n, unsigned const * const iA, unsigned const * const jA,
	double const * const A, double const * const x, double* y,
	unsigned num_of_pthreads);

void amuxCRSParallelPThreads (double a,
	unsigned n, unsigned const * const iA, unsigned const * const jA,
	double const * const A, double const * const x, double* y,
	unsigned num_of_pthreads, unsigned const*const workload_intervals);

#ifdef _OPENMP
template<typename FP_TYPE, typename IDX_TYPE>
void amuxCRSParallelOpenMP (FP_TYPE a, unsigned n,
    IDX_TYPE const * const __restrict__ iA,
    IDX_TYPE const * const __restrict__ jA, FP_TYPE const * const A,
    FP_TYPE const * const __restrict__ x, FP_TYPE* __restrict__ y)
{
	OPENMP_LOOP_TYPE i;
    IDX_TYPE j;
    FP_TYPE t;
	{
#pragma omp parallel for private(i, j, t)
		for (i = 0; i < n; i++) {
			const IDX_TYPE end(iA[i + 1]);
			t = A[iA[i]] * x[jA[iA[i]]];
			for (j = iA[i]+1; j < end; j++) {
				t += A[j] * x[jA[j]];
			}
            y[i] = t * a;
		}
	}
}
#endif

void amuxCRSSym (double a,
	unsigned n, unsigned const * const iA, unsigned const * const jA,
        double const * const A, double const * const x, double* y);

} // end namespace MathLib

#endif
