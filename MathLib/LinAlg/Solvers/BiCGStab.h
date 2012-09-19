/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file BiCGStab.h
 *
 * Created on 2011-10-04 by Thomas Fischer
 */

#ifndef BICGSTAB_H_
#define BICGSTAB_H_

#include "blas.h"
#include "../Sparse/CRSMatrix.h"
#include "../Sparse/CRSMatrixDiagPrecond.h"

namespace MathLib {

unsigned BiCGStab(CRSMatrix<double, unsigned> const& A, double* const b, double* const x,
                  double& eps, unsigned& nsteps);

} // end namespace MathLib

#endif /* BICGSTAB_H_ */
