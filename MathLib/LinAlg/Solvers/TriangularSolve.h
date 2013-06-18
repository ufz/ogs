/**
 * \file
 * \author Thomas Fischer
 * \date   2011-05-06
 * \brief  Definition of triangular solver functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TRIANGULARSOLVE_H_
#define TRIANGULARSOLVE_H_

#include "../Dense/DenseMatrix.h"

namespace MathLib {

/**
 * solves the \f$n \times n\f$ triangular linear system \f$L \cdot y = b\f$,
 * assumes \f$L_{ii} = 1.0\f$, \f$i=1,...,n\f$, \f$b\f$ is destroyed
 * @param L the lower triangular matrix
 * @param b at beginning the right hand side vector, at the end the solution vector
 */
template <typename FP_T>
void forwardSolve (const DenseMatrix <FP_T> &L, FP_T* b);

/**
 * solves the \f$n \times n\f$ triangular linear system \f$U \cdot x=y\f$,
 * \f$U\f$, where \f$U\f$ is a upper triangular matrix.
 * @param U upper triangular matrix
 * @param y at beginning the right hand side, at the end the solution
 */
template <typename FP_T>
void backwardSolve (const DenseMatrix <FP_T> &U, FP_T* y);

// backwardSolve mat * x = y, mat ... upper triangular matrix
/**
 * backward solve the system of linear equations \f$ U \cdot x = y\f$,
 * where \f$U\f$ is a upper triangular matrix
 * @param mat the upper triangular matrix
 * @param x the solution of the system of linear equations
 * @param b the right hand side
 */
template <typename FP_T>
void backwardSolve ( DenseMatrix<FP_T> const& mat, FP_T* x, FP_T* b);

} // end namespace MathLib

#include "TriangularSolve.tpp"

#endif /* TRIANGULARSOLVE_H_ */
