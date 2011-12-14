/*
 * CGParallel.cpp
 *
 *  Created on: Dec 2, 2011
 *      Author: TF
 */

#include <limits>

#include <omp.h>

#include "MathTools.h"
#include "blas.h"
#include "../Sparse/CRSMatrix.h"
#include "../Sparse/CRSMatrixDiagPrecond.h"

// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//      x  --  approximate solution to Ax = b
// nsteps  --  the number of iterations performed before the
//             tolerance was reached
//    eps  --  the residual after the final iteration

namespace MathLib {

#ifdef _OPENMP
unsigned CGParallel(CRSMatrix<double,unsigned> const * mat, double const * const b,
		double* const x, double& eps, unsigned& nsteps, unsigned num_threads)
{
	omp_set_num_threads(num_threads);
	unsigned N = mat->getNRows();
	double * __restrict__ p(new double[N]);
	double * __restrict__ q(new double[N]);
	double * __restrict__ r(new double[N]);
	double * __restrict__ rhat(new double[N]);
	double rho, rho1 = 0.0;

//	p = new double[4* N];
//	q = p + N;
//	r = q + N;
//	rhat = r + N;

	double nrmb = sqrt(scpr(b, b, N));
	if (nrmb < std::numeric_limits<double>::epsilon()) {
		blas::setzero(N, x);
		eps = 0.0;
		nsteps = 0;
		delete[] p;
		return 0;
	}

	// r0 = b - Ax0
	mat->amux(D_MONE, x, r);
	for (unsigned k(0); k < N; k++) {
		r[k] = b[k] - r[k];
	}

	double resid = blas::nrm2(N, r);
	if (resid <= eps * nrmb) {
		eps = resid / nrmb;
		nsteps = 0;
		delete[] p;
		return 0;
	}

	unsigned k;
	for (unsigned l = 1; l <= nsteps; ++l) {
#ifndef NDEBUG
		std::cout << "Step " << l << ", resid=" << resid / nrmb << std::endl;
#endif
		// r^ = C r
		// rhat = r
//		blas::copy(N, r, rhat);
		#pragma omp parallel for
		for (k = 0; k < N; k++) {
			rhat[k] = r[k];
		}
		mat->precondApply(rhat);

		// rho = r * r^;
		rho = scpr(r, rhat, N, num_threads);

		if (l > 1) {
			double beta = rho / rho1;
			// p = r^ + beta * p
			#pragma omp parallel for
			for (k = 0; k < N; k++) {
				p[k] = rhat[k] + beta * p[k];
			}
		} else {
//			blas::copy(N, rhat, p);
			#pragma omp parallel for
			for (k = 0; k < N; k++) {
				p[k] = rhat[k];
			}
		}

		// q = Ap
//		blas::setzero(N, q);
		#pragma omp parallel for
		for (k = 0; k < N; k++) {
			q[k] = 0.0;
		}
		mat->amux(D_ONE, p, q);

		// alpha = rho / p*q
		double alpha = rho / scpr(p, q, N, num_threads);

		// x += alpha * p
		#pragma omp parallel for
		for (k = 0; k < N; k++) {
			x[k] += alpha * p[k];
		}

		// r -= alpha * q
		#pragma omp parallel for
		for (k = 0; k < N; k++) {
			r[k] -= alpha * q[k];
		}

		resid = sqrt(scpr(r, r, N, num_threads));

		if (resid <= eps * nrmb) {
			eps = resid / nrmb;
			nsteps = l;
			delete[] p;
			return 0;
		}

		rho1 = rho;
	}
	eps = resid / nrmb;
	delete[] p;
	return 1;
}
#endif

} // end of namespace MathLib
