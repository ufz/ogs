/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of Lis utility functions.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISTOOLS_H_
#define LISTOOLS_H_

#include <vector>

namespace MathLib
{
class LisMatrix;
class LisVector;

/**
 * apply known solutions to a system of linear equations
 *
 * This function introduces the constants into the system by the penalty method.
 *
 * @param A                 Coefficient matrix
 * @param b                 RHS vector
 * @param _vec_knownX_id    a vector of known solution entry IDs
 * @param _vec_knownX_x     a vector of known solutions
 * @param penalty_scaling value for scaling some matrix and right hand side
 * entries to enforce some conditions
 */
void applyKnownSolution(LisMatrix &A, LisVector &b, const std::vector<std::size_t> &_vec_knownX_id,
		const std::vector<double> &_vec_knownX_x, double penalty_scaling = 1e+10);

} // MathLib

#endif //LISTOOLS_H_

