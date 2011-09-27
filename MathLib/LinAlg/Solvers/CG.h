/*
 * CG.h
 *
 *  Created on: Sep 27, 2011
 *      Author: TF
 */

#include "blas.h"
#include "CG.h"
#include <limits>
#include "CRSMatrix.h"
#include "CRSMatrixDiagPrecond.h"
#include "sparse.h"

#include <fstream>

#ifndef NDEBUG
#include <iostream>
#endif

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


unsigned CG(CRSMatrix<double> const * mat, double const * const b,
	double* const x, double& eps, unsigned& nsteps, unsigned num_threads)
{
  unsigned N = mat->getNRows();
  double *p, *q, *r, *rhat, rho, rho1 = 0.0;

  p = new double[4*N];
  q = p + N;
  r = q + N;
  rhat = r + N;

  double nrmb = sqrt(scpr(N, b, b, num_threads));
  if (nrmb < std::numeric_limits<double>::epsilon()) {
    blas::setzero(N, x);
    eps = 0.0;
    nsteps = 0;
    delete [] p;
    return 0;
  }

  // r0 = b - Ax0
  mat->amux (D_MONE, x, r);
  for (unsigned k(0); k<N; k++) {
  	r[k] = b[k] - r[k];
  }

  std::ofstream out ("residuen.txt");

  double resid = blas::nrm2(N, r);
  out << "0\t" << resid/nrmb << std::endl;
  if (resid <= eps*nrmb) {
    eps = resid / nrmb;
    nsteps = 0;
    delete [] p;
    out.close();
    return 0;
  }

  for (unsigned l=1; l<=nsteps; ++l) {

	  out << l << "\t" << resid/nrmb << std::endl;
#ifndef NDEBUG
    std::cout << "Step " << l << ", resid=" << resid/nrmb << std::endl;
#endif

    // r^ = C r
    blas::copy(N, r, rhat);
    mat->precondApply(rhat);

    // rho = r * r^;
    rho = scpr(N, r, rhat, num_threads);

    if (l>1) {
	double beta = rho / rho1;
	// p = r^ + beta * p
	unsigned k;
	omp_set_num_threads (num_threads);
	#pragma omp parallel for
	for(k=0; k<N; k++) {
		p[k] = rhat[k] + beta * p[k];
	}
    }
    else
      blas::copy(N, rhat, p);

    // q = Ap
    blas::setzero(N, q);
    mat->amux (D_ONE, p, q);

    // alpha = rho / p*q
    double alpha = rho / scpr(N, p, q);

    // x += alpha * p
    blas::axpy(N, alpha, p, x);

    // r -= alpha * q
    blas::axpy(N, -alpha, q, r);

    resid = sqrt(scpr(N,r,r));

    if (resid<=eps*nrmb) {
      eps = resid/nrmb;
      nsteps = l;
      delete [] p;
      return 0;
    }

    rho1 = rho;
  }
  out.close();
  eps = resid / nrmb;
  delete [] p;
  return 1;
}


