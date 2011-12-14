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

#ifdef _OPENMP
template<typename FP_TYPE, typename IDX_TYPE>
void amuxCRSParallelOpenMP (FP_TYPE a,
				unsigned n, IDX_TYPE const * const __restrict__ iA, IDX_TYPE const * const __restrict__ jA,
				FP_TYPE const * const A, FP_TYPE const * const __restrict__ x, FP_TYPE* __restrict__ y)
{
	unsigned i;
	{
#pragma omp parallel for
		for (i = 0; i < n; i++) {
			const IDX_TYPE end(iA[i + 1]);
			y[i] = A[iA[i]] * x[jA[iA[i]]];
			for (IDX_TYPE j(iA[i]+1); j < end; j++) {
				y[i] += A[j] * x[jA[j]];
			}
			y[i] *= a;
		}
	}
}
#endif

void amuxCRSSym (double a,
	unsigned n, unsigned const * const iA, unsigned const * const jA,
        double const * const A, double const * const x, double* y);

} // end namespace MathLib

#endif
