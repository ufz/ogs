/**
 * \file
 * \author Thomas Fischer
 * \date   2012-01-20
 * \brief  Definition of the CRSMatrixReorderedOpenMP class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
