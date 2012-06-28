/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file GMRes.cpp
 *
 * Created on 2011-10-04 by Thomas Fischer
 */

#include "GMRes.h"

#include <cmath>
#include <limits>
#include "blas.h"

namespace MathLib {

static void genPlRot(double dx, double dy, double& cs, double& sn)
{
	if (dy <= std::numeric_limits<double>::epsilon()) {
		cs = 1.0;
		sn = 0.0;
	} else if (fabs(dy) > fabs(dx)) {
		const double tmp = dx / dy;
		sn = 1.0 / sqrt(1.0 + tmp * tmp);
		cs = tmp * sn;
	} else {
		const double tmp = dy / dx;
		cs = 1.0 / sqrt(1.0 + tmp * tmp);
		sn = tmp * cs;
	}
}

inline void applPlRot(double& dx, double& dy, double cs, double sn)
{
	const double tmp = cs * dx + sn * dy;
	dy = cs * dy - sn * dx;
	dx = tmp;
}

// solve H y = s and update x += MVy
static void update(const CRSMatrix<double,unsigned>& A, unsigned k, double* H,
		unsigned ldH, double* s, double* V, double* x)
{
	const size_t n(A.getNRows());
	double *y = new double[k];
	double *xh = new double[n];
	blas::copy(k, s, y);
	int inf;

	dtrtrs_(JOB_STR + 5, JOB_STR, JOB_STR, &k, &N_ONE, H, &ldH, y, &k, &inf);
	assert(inf == 0);

	// x += M V y
	blas::setzero(n, xh);
	blas::gemva(n, k, D_ONE, V, y, xh);
	A.precondApply(xh);
	blas::add(n, xh, x);

	delete[] xh;
	delete[] y;
}

unsigned GMRes(const CRSMatrix<double,unsigned>& A, double* const b, double* const x,
		double& eps, unsigned m, unsigned& nsteps)
{
	double resid;
	unsigned j = 1;

	const size_t n (A.getNRows());

	double *r = new double[2*n + (n + m + 4) * (m + 1)]; // n
	double *V = r + n; // n x (m+1)
	double *H = V + n * (m + 1); // m+1 x m
	double *cs = H + (m + 1) * m; // m+1
	double *sn = cs + m + 1; // m+1
	double *s = sn + m + 1; // m+1
	double *xh = s + m + 1; // m+1

	// normb = norm(b)
	double normb = blas::nrm2(n, b);
	if (normb == 0.0) {
		blas::setzero(n, x);
		eps = 0.0;
		nsteps = 0;
		delete[] r;
		return 0;
	}

	// r = b - Ax
	blas::copy(n, b, r);
	A.amux(D_MONE, x, r);

	double beta = blas::nrm2(n, r);

	if ((resid = beta / normb) <= eps) {
		eps = resid;
		nsteps = 0;
		delete[] r;
		return 0;
	}

	while (j <= nsteps) {
		blas::copy(n, r, V); // v0 first orthonormal vector
		blas::scal(n, 1.0 / beta, V);

		s[0] = beta;
		blas::setzero(m, s + 1);

		for (unsigned i = 0; i < m && j <= nsteps; i++, j++) {

			// w = A M * v[i];
			blas::copy(n, V + i * n, xh);
			A.precondApply(xh);
			blas::setzero(n, V + (i + 1) * n);
			A.amux(D_ONE, xh, V + (i + 1) * n);

			for (unsigned k = 0; k <= i; k++) {
				H[k + i * (m + 1)] = blas::scpr(n, V + (i + 1) * n, V + k * n);
				blas::axpy(n, -H[k + i * (m + 1)], V + k * n, V + (i + 1) * n);
			}

			H[i * (m + 2) + 1] = blas::nrm2(n, V + (i + 1) * n);
			blas::scal(n, 1.0 / H[i * (m + 2) + 1], V + (i + 1) * n);

			// apply old Givens rotations to the last column in H
			for (unsigned k = 0; k < i; k++)
				applPlRot(H[k + i * (m + 1)], H[k + 1 + i * (m + 1)], cs[k],
						sn[k]);

			// generate new Givens rotation which eleminates H[i*(m+2)+1]
			genPlRot(H[i * (m + 2)], H[i * (m + 2) + 1], cs[i], sn[i]);
			// apply it to H and s
			applPlRot(H[i * (m + 2)], H[i * (m + 2) + 1], cs[i], sn[i]);
			applPlRot(s[i], s[i + 1], cs[i], sn[i]);

			if ((resid = fabs(s[i + 1] / normb)) < eps) {
				update(A, i + 1, H, m + 1, s, V, x);
				eps = resid;
				nsteps = j;
				delete[] r;
				return 0;
			}
#ifndef NDEBUG
			std::cout << "Step " << j << ", resid=" << resid << std::endl;
#endif
		}

		update(A, m, H, m + 1, s, V, x);

		// r = b - A x;
		blas::copy(n, b, r);
		A.amux(D_MONE, x, r);
		beta = blas::nrm2(n, r);

		if ((resid = beta / normb) < eps) {
			eps = resid;
			nsteps = j;
			delete[] r;
			return 0;
		}
	}

	eps = resid;
	delete[] r;
	return 1;
}

} // end namespace MathLib
