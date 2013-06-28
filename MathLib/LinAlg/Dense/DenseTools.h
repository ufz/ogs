/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Declaration of dense matrix and vector utility functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DENSETOOLS_H_
#define DENSETOOLS_H_

#include <vector>
#include "DenseMatrix.h"
#include "DenseVector.h"

namespace MathLib
{

/**
 * Apply known solutions to a system of linear equations.
 *
 * This function introduces the given constrains by diagonalizing a coefficient matrix.
 * Symmetry of the matrix is preserved.
 *
 * @param A                Coefficient matrix
 * @param b                RHS vector
 * @param vec_knownX_id    a vector of known solution entry IDs
 * @param vec_knownX_x     a vector of known solutions
 */
void applyKnownSolution(DenseMatrix<double> &A, DenseVector<double> &b,
		const std::vector<std::size_t> &vec_knownX_id, const std::vector<double> &vec_knownX_x);

/**
 * Apply known solutions to a system of linear equations \f$A x = b\f$.
 *
 * This function introduces the given constrain into the system of linear
 * equations. For this purpose it modifies the entries of \f$A\f$ in the
 * \f$k\f$-th row and column, i.e. all those entries are set to zero except
 * the diagonal entry that is set to one. The right
 * hand side \f$b\f$ is modified, too. The symmetry of \f$A\f$ is preserved.
 *
 * @param A         Coefficient matrix
 * @param b         RHS vector
 * @param row_id    a known solution entry ID
 * @param val       a known solution
 */
void applyKnownSolution(DenseMatrix<double> &A, DenseVector<double> &b, std::size_t row_id,
		double val);

} // MathLib

#endif //DENSETOOLS_H_

