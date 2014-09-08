/**
 * \file
 * \author Thomas Fischer
 * \date   2012-01-12
 * \brief  Implementation of the CRSMatrixReorderedOpenMP class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifdef _OPENMP
#include "LinAlg/Sparse/NestedDissectionPermutation/CRSMatrixReorderedOpenMP.h"
#include "LinAlg/Sparse/amuxCRS.h"

namespace MathLib {

CRSMatrixReorderedOpenMP::CRSMatrixReorderedOpenMP(unsigned n, unsigned *iA, unsigned *jA, double* A) :
	CRSMatrixReordered(n, iA, jA, A)
{
}

CRSMatrixReorderedOpenMP::~CRSMatrixReorderedOpenMP()
{}

void CRSMatrixReorderedOpenMP::amux(double d, double const * const x, double *y) const
{
	amuxCRSParallelOpenMP(d, this->_n_rows, this->_row_ptr, this->_col_idx, this->_data, x, y);
}

} // end namespace MathLib

#endif // _OPENMP
