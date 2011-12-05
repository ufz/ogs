#ifndef AMUXCRS_H
#define AMUXCRS_H

namespace MathLib {

template<typename FP_TYPE, typename IDX_TYPE>
void amuxCRS(FP_TYPE a, IDX_TYPE n, IDX_TYPE const * const iA, IDX_TYPE const * const jA,
				FP_TYPE const * const A, FP_TYPE const * const x, FP_TYPE* y)
{
	for (IDX_TYPE i(0); i < n; i++) {
		y[i] = 0.0;
		const IDX_TYPE end(iA[i + 1]);
		for (IDX_TYPE j(iA[i]); j < end; j++) {
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
				unsigned n, IDX_TYPE const * const iA, IDX_TYPE const * const jA,
				FP_TYPE const * const A, FP_TYPE const * const x, FP_TYPE* y,
				unsigned num_of_omp_threads)
{
	unsigned i;
	{
#pragma omp parallel for
		for (i = 0; i < n; i++) {
			FP_TYPE temp (0.0);
			const IDX_TYPE end(iA[i + 1]);
			for (IDX_TYPE j(iA[i]); j < end; j++) {
				temp += A[j] * x[jA[j]];
			}
			y[i] = a*temp;
		}
	}
}
#endif

void amuxCRSSym (double a,
	unsigned n, unsigned const * const iA, unsigned const * const jA,
        double const * const A, double const * const x, double* y);

} // end namespace MathLib

#endif
