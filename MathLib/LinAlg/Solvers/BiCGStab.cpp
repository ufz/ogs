/**
 * \file
 * \author Thomas Fischer
 * \date   2011-10-04
 * \brief  Implementation of the BiCGStab function.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BiCGStab.h"

#include "MathTools.h"
#include "blas.h"

namespace MathLib {

unsigned BiCGStab(CRSMatrix<double, unsigned> const& A, double* const b, double* const x,
		double& eps, unsigned& nsteps)
{
	const unsigned N(A.getNRows());
	double *v (new double[8* N]);
	double *p (v + N);
	double *phat (p + N);
	double *s (phat + N);
	double *shat (s + N);
	double *t (shat + N);
	double *r (t + N);
	double *r0 (r + N);
	double resid;

	// normb = |b|
	double nrmb = blas::nrm2(N, b);
	if (nrmb < D_PREC) nrmb = D_ONE;

	// r = r0 = b - A x0
	blas::copy(N, b, r0);
	A.amux(D_MONE, x, r0);
	blas::copy(N, r0, r);

	resid = blas::nrm2(N, r) / nrmb;

	if (resid < eps) {
		eps = resid;
		nsteps = 0;
		delete[] v;
		return 0;
	}

	double alpha = D_ZERO, omega = D_ZERO, rho2 = D_ZERO;

	for (unsigned l = 1; l <= nsteps; ++l) {
		// rho1 = r0 * r
		const double rho1 = blas::scpr(N, r0, r);
		if (fabs(rho1) < D_PREC) {
			eps = blas::nrm2(N, r) / nrmb;
			delete[] v;
			return 2;
		}

		if (l == 1)
			blas::copy(N, r, p); // p = r
		else {
//			blas::axpy(N, -omega, v, p); // p = (p-omega v)*beta+r
			const double beta = rho1 * alpha / (rho2 * omega);
//			blas::scal(N, beta, p);
//			blas::axpy(N, D_ONE, r, p);
			// p = (p-omega v)*beta+r
			for (unsigned k(0); k<N; k++) {
				p[k] = (p[k] - omega * v[k]) * beta + r[k];
			}
		}

		// p^ = C p
		blas::copy(N, p, phat);
		A.precondApply(phat);
		// v = A p^
		blas::setzero(N, v);
		A.amux(D_ONE, phat, v);

		alpha = rho1 / blas::scpr(N, r0, v);

		// s = r - alpha v
//		blas::copy(N, r, s);
//		blas::axpy(N, -alpha, v, s);
		for (unsigned k(0); k<N; k++) {
			s[k] = r[k] - alpha * v[k];
		}

		resid = blas::nrm2(N, s) / nrmb;
#ifndef NDEBUG
		std::cout << "Step " << l << ", resid=" << resid << std::endl;
#endif
		if (resid < eps) {
			// x += alpha p^
			blas::axpy(N, alpha, phat, x);
			eps = resid;
			nsteps = l;
			delete[] v;
			return 0;
		}

		// s^ = C s
		blas::copy(N, s, shat);
		A.precondApply(shat);

		// t = A s^
		blas::setzero(N, t);
		A.amux(D_ONE, shat, t);

		// omega = t*s / t*t
		omega = blas::scpr(N, t, s) / blas::scpr(N, t, t);

		// x += alpha p^ + omega s^
		blas::axpy(N, alpha, phat, x);
		blas::axpy(N, omega, shat, x);

		// r = s - omega t
		blas::copy(N, s, r);
		blas::axpy(N, -omega, t, r);

		rho2 = rho1;

		resid = blas::nrm2(N, r) / nrmb;

		if (resid < eps) {
			eps = resid;
			nsteps = l;
			delete[] v;
			return 0;
		}

		if (fabs(omega) < D_PREC) {
			eps = resid;
			delete[] v;
			return 3;
		}
	}

	eps = resid;
	delete[] v;
	return 1;
}

} // end namespace MathLib
