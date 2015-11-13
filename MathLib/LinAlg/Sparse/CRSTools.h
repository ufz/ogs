/**
 * \brief Integration of Dirichlet boundary conditions into the system of
 * linear equations.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CRSTOOLS_H_
#define CRSTOOLS_H_

#include <vector>

namespace MathLib
{

template<typename FP_TYPE, typename IDX_TYPE> class CRSMatrix;
template<typename IDX_TYPE> class DenseVector;

/**
 * Integrate Dirichlet boundary conditions into a system of linear equations.
 *
 * This function introduces the constants into the system by setting
 * appropriated row and column entries of the matrix to zero (except the
 * diagonal entries) and modifying values within the right hand side vector.
 *
 * @param mat Coefficient matrix
 * @param rhs RHS vector
 * @param rows a vector of known solution entry IDs
 * @param vals a vector of known solutions
 */
template <typename VEC_T, typename FP_TYPE = double>
void applyKnownSolution(CRSMatrix<FP_TYPE, typename VEC_T::IndexType>*& mat,
	VEC_T &rhs, std::vector<std::size_t> const& rows,
	std::vector<FP_TYPE> const& vals);

} // MathLib

#include "CRSTools-impl.h"

#endif // CRSTOOLS_H_

