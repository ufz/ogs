/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CG.cpp
 *
 * Created on 2011-09-27 by Thomas Fischer
 */

#include <limits>

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

unsigned CG(CRSMatrix<double,unsigned> const * mat, double const * const b,
		double* const x, double& eps, unsigned& nsteps)
{
	unsigned N = mat->getNRows();
	double *p, *q, *r, *rhat, rho, rho1 = 0.0;

	p = new double[4* N];
	q = p + N;
	r = q + N;
	rhat = r + N;

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

	for (unsigned l = 1; l <= nsteps; ++l) {
#ifndef NDEBUG
		std::cout << "Step " << l << ", resid=" << resid / nrmb << std::endl;
#endif
		// r^ = C r
		blas::copy(N, r, rhat);
		mat->precondApply(rhat);

		// rho = r * r^;
		rho = scpr(r, rhat, N); // num_threads);

		if (l > 1) {
			double beta = rho / rho1;
			// p = r^ + beta * p
			unsigned k;
			for (k = 0; k < N; k++) {
				p[k] = rhat[k] + beta * p[k];
			}
		} else blas::copy(N, rhat, p);

		// q = Ap
		blas::setzero(N, q);
		mat->amux(D_ONE, p, q);

		// alpha = rho / p*q
		double alpha = rho / scpr(p, q, N);

		// x += alpha * p
		blas::axpy(N, alpha, p, x);

		// r -= alpha * q
		blas::axpy(N, -alpha, q, r);

		resid = sqrt(scpr(r, r, N));

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

} // end namespace MathLib
