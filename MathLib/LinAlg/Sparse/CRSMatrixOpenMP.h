/**
 * \file
 * \author Thomas Fischer
 * \date   2011-08-08
 * \brief  Definition of the CRSMatrixOpenMP class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CRSMATRIXOPENMP_H_
#define CRSMATRIXOPENMP_H_

#ifdef _OPENMP
#include <string>

#include "CRSMatrix.h"
#include "amuxCRS.h"

namespace MathLib {

template<typename FP_TYPE, typename IDX_TYPE> class CRSMatrixOpenMP : public CRSMatrix<FP_TYPE, IDX_TYPE> {
public:
	CRSMatrixOpenMP(std::string const &fname) :
			CRSMatrix<FP_TYPE, IDX_TYPE>(fname)
	{}

	CRSMatrixOpenMP(unsigned n, IDX_TYPE *iA, IDX_TYPE *jA, FP_TYPE* A) :
		CRSMatrix<FP_TYPE, IDX_TYPE>(n, iA, jA, A)
	{}

	CRSMatrixOpenMP(unsigned n1) :
		CRSMatrix<FP_TYPE, IDX_TYPE>(n1)
	{}

	virtual ~CRSMatrixOpenMP()
	{}

	virtual void amux(FP_TYPE const d, FP_TYPE const * const x, FP_TYPE *y) const
	{
		amuxCRSParallelOpenMP(d, this->_n_rows, CRSMatrix<FP_TYPE,IDX_TYPE>::_row_ptr, CRSMatrix<FP_TYPE,IDX_TYPE>::_col_idx, CRSMatrix<FP_TYPE,IDX_TYPE>::_data, x, y);
	}
};

} // end namespace MathLib

#endif // _OPENMP

#endif /* CRSMATRIXOPENMP_H_ */
