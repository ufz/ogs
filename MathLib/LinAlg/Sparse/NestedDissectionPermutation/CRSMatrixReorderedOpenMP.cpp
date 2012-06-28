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
 * \file CRSMatrixReorderedOpenMP.cpp
 *
 * Created on 2012-01-12 by Thomas Fischer
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
	amuxCRSParallelOpenMP(d, MatrixBase::_n_rows, CRSMatrix<double,unsigned>::_row_ptr, CRSMatrix<double,unsigned>::_col_idx, CRSMatrix<double,unsigned>::_data, x, y);
}

} // end namespace MathLib

#endif // _OPENMP
