/*
 * GMRes.h
 *
 *  Created on: Oct 4, 2011
 *      Author: TF
 */

#ifndef GMRES_H_
#define GMRES_H_

#include "../Sparse/CRSMatrix.h"

namespace MathLib {

unsigned GMRes(const CRSMatrix<double,unsigned>& mat, double* const b, double* const x,
                        double& eps, unsigned m, unsigned& steps);

} // end namespace MathLib

#endif /* GMRES_H_ */
