/*
 * CRSMatrixReordered.h
 *
 *  Created on: Jan 3, 2012
 *      Author: TF
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
