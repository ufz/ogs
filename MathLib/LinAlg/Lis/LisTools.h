/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of Lis utility functions.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISTOOLS_H_
#define LISTOOLS_H_

#include <vector>

#include "LisMatrix.h" // for LisMatrix::IndexType

namespace MathLib
{
class LisVector;

/**
 * Integrate Dirichlet boundary conditions into a system of linear equations.
 *
 * This function introduces the constants into the system by setting
 * appropriated row and column entries of the matrix to zero (except the
 * diagonal entries) and modifying values within the right hand side vector.
 *
 * @param eqsA                 Coefficient matrix
 * @param eqsRHS                 RHS vector
 * @param rows a vector of known solution entry IDs
 * @param vals a vector of known solutions
 */
void applyKnownSolution(LisMatrix &eqsA, LisVector &eqsRHS, LisVector &/*eqsX*/,
	const std::vector<LisMatrix::IndexType> &rows,
	const std::vector<double> &vals);

} // MathLib

#endif //LISTOOLS_H_

