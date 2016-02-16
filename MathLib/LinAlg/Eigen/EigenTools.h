/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENTOOLS_H_
#define EIGENTOOLS_H_

#include <vector>

#include "EigenMatrix.h" // for EigenMatrix::IndexType

namespace MathLib
{
class EigenVector;

/**
 * apply known solutions to a system of linear equations
 *
 * @param A                 Coefficient matrix
 * @param b                 RHS vector
 * @param _vec_knownX_id    a vector of known solution entry IDs
 * @param _vec_knownX_x     a vector of known solutions
 * @param penalty_scaling value for scaling some matrix and right hand side
 * entries to enforce some conditions
 */
void applyKnownSolution(EigenMatrix &A, EigenVector &b, EigenVector &/*x*/,
		const std::vector<EigenMatrix::IndexType> &_vec_knownX_id,
		const std::vector<double> &_vec_knownX_x, double penalty_scaling = 1e+10);

inline
void applyKnownSolution(Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::VectorXd &/*x*/,
		const std::vector<Eigen::MatrixXd::Index> &_vec_knownX_id,
		const std::vector<double> &_vec_knownX_x, double penalty_scaling = 1e+10)
{
	// TODO implement
	(void) A; (void) b; (void) _vec_knownX_id; (void) _vec_knownX_x;
	(void) penalty_scaling;
}

} // MathLib

#endif //EIGENTOOLS_H_

