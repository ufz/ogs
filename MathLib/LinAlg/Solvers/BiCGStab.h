/*
 * BiCGStab.h
 *
 *  Created on: Oct 4, 2011
 *      Author: TF
 */

#ifndef BICGSTAB_H_
#define BICGSTAB_H_

#include "blas.h"
#include "../Sparse/CRSMatrix.h"
#include "../Sparse/CRSMatrixDiagPrecond.h"

namespace MathLib {

unsigned BiCGStab(CRSMatrix<double> const& A, double* const b, double* const x,
                  double& eps, unsigned& nsteps);

} // end namespace MathLib

#endif /* BICGSTAB_H_ */
