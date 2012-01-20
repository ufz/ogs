/*
 * CRSMatrixReorderedOpenMP.h
 *
 *  Created on: Jan 20, 2012
 *      Author: TF
 */

#ifndef CRSMATRIXREORDEREDOPENMP_H_
#define CRSMATRIXREORDEREDOPENMP_H_

#ifdef _OPENMP
#include "LinAlg/Sparse/NestedDissectionPermutation/CRSMatrixReordered.h"

namespace MathLib {

class CRSMatrixReorderedOpenMP : public CRSMatrixReordered {
public:
	CRSMatrixReorderedOpenMP(unsigned n, unsigned *iA, unsigned *jA, double* A);
	virtual ~CRSMatrixReorderedOpenMP();

	virtual void amux(double d, double const * const x, double *y) const;
};

}
#endif // _OPENMP

#endif /* CRSMATRIXREORDEREDOPENMP_H_ */
