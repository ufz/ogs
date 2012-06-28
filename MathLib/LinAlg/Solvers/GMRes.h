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
 * \file GMRes.h
 *
 * Created on 2011-10-04 by Thomas Fischer
 */

#ifndef GMRES_H_
#define GMRES_H_

#include "../Sparse/CRSMatrix.h"

namespace MathLib {

unsigned GMRes(const CRSMatrix<double,unsigned>& mat, double* const b, double* const x,
                        double& eps, unsigned m, unsigned& steps);

} // end namespace MathLib

#endif /* GMRES_H_ */
