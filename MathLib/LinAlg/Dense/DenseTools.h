/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of dense matrix and vector utility functions.
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
#include "Matrix.h"
#include "Vector.h"

namespace MathLib
{

/**
 * apply known solutions to a system of linear equations
 *
 * This function introduces the given constrain by diagonalizing a coefficient matrix.
 * Symmetricity of the matrix is preserved.
 *
 * @param A                 Coefficient matrix
 * @param b                 RHS vector
 * @param _vec_knownX_id    a vector of known solution entry IDs
 * @param _vec_knownX_x     a vector of known solutions
 */
void applyKnownSolution(Matrix<double> &A, Vector<double> &b, const std::vector<std::size_t> &_vec_knownX_id, const std::vector<double> &_vec_knownX_x);

/**
 * apply known solutions to a system of linear equations
 *
 * This function introduces the given constrain by diagonalizing a coefficient matrix.
 * Symmetricity of the matrix is preserved.
 *
 * @param A         Coefficient matrix
 * @param b         RHS vector
 * @param row_id    a known solution entry ID
 * @param x         a known solution
 */
void applyKnownSolution(MathLib::Matrix<double> &A, MathLib::Vector<double> &b, size_t row_id, double x);

} // MathLib

#endif //DENSETOOLS_H_

