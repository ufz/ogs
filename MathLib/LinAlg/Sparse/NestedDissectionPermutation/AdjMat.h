/**
 * \file
 * \author Thomas Fischer
 * \date   2012-01-02
 * \brief  Definition of the AdjMat class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ADJMAT_H_
#define ADJMAT_H_

#include "LinAlg/Sparse/CRSMatrix.h"

namespace MathLib {

/** AdjMat implements an adjacency matrix used to represent the (possible
   weighted) edges of a graph. Due to adjacency matrix are mostly sparse
   the implementation makes use of compressed row storage format.
*/
class AdjMat: public CRSMatrix<unsigned, unsigned>
{
public:
	/** constructor with data elements in A
	 * @param s size of the quadratic matrix
	 * @param iA array of length s+1 which holds pointers in jA,
	 * iA[k] points to the first non-zero column-entry of row k,
	 * iA[k]-1 points accordingly to the last non-zero column-entry of row k,
	 * the last entry of iA (iA[s]) takes the number of non zero entries(nnz)
	 * @param jA array of length nnz, each entry is a colum-index
	 * @param A  data-array of length nnz of type unsigned (weights)
	 */
	AdjMat(unsigned s, unsigned *iA, unsigned *jA, unsigned *A = NULL) :
		CRSMatrix<unsigned, unsigned> (s, iA, jA, A)
	{}

	/**
	 * destructor
	 */
	virtual ~AdjMat()
	{}

	/**
	 *
	 */
	void makeSymmetric();

	/** getMat returns the (possibly reducible) block [beg,end-1] x [beg,end-1]
	 * respecting the permutation.
	 * @param beg index of first row/column, it is supposed that 0 <= beg <= n
	 * @param end index one after last row/column, it is supposed that beg <= end <= n
	 * @param op_perm permutation -> original
	 * @param po_perm original -> permutation
	 * @return pointer to an AdjMat object
	 */
	AdjMat* getMat(unsigned beg, unsigned end, unsigned const* const op_perm,
			unsigned const* const po_perm) const;

private:
	DISALLOW_COPY_AND_ASSIGN(AdjMat);
};

}

#endif /* ADJMAT_H_ */
