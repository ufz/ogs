/**
 * \file
 * \author Thomas Fischer
 * \date   2011-10-04
 * \brief  Definition of the BiCGStab function.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
