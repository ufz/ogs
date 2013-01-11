/**
 * \file
 * \author Thomas Fischer
 * \date   2011-10-04
 * \brief  Definition of the GMRes function.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GMRES_H_
#define GMRES_H_

#include "../Sparse/CRSMatrix.h"

namespace MathLib {

unsigned GMRes(const CRSMatrix<double,unsigned>& mat, double* const b, double* const x,
                        double& eps, unsigned m, unsigned& steps);

} // end namespace MathLib

#endif /* GMRES_H_ */
