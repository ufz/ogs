/**
 * \file
 * \author Thomas Fischer
 * \date   2012-01-03
 * \brief  Definition of the CRSMatrixReordered class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CRSMATRIXREORDERED_H_
#define CRSMATRIXREORDERED_H_

#include "LinAlg/Sparse/CRSMatrix.h"

namespace MathLib {

class CRSMatrixReordered: public MathLib::CRSMatrix<double,unsigned>
{
public:
	CRSMatrixReordered(std::string const &fname);
	CRSMatrixReordered(unsigned n, unsigned *iA, unsigned *jA, double* A);
	virtual ~CRSMatrixReordered();
	void reorderMatrix(unsigned const*const op_perm, unsigned const*const po_perm);
};

}

#endif /* CRSMATRIXREORDERED_H_ */
